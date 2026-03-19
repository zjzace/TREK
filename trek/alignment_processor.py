#!/usr/bin/env python3
"""
Alignment Processor: Handle minimap2 alignment and read-to-transcript assignment
"""

import re
import subprocess
import logging
import pysam
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from dataclasses import dataclass
from tqdm import tqdm

logger = logging.getLogger(__name__)


@dataclass
class AlignedRead:
    """Represents an aligned read with assignment info"""
    read_id: str
    chromosome: str
    start: int  # 0-based
    end: int    # 0-based
    strand: str
    cigar: str
    mapq: int
    
    def has_large_clips(self, max_clip_size: int = 20) -> bool:
        """
        Check if read has large soft/hard clips at either end
        
        Args:
            max_clip_size: Maximum allowed clip size
            
        Returns:
            True if clips exceed max_clip_size
        """
        cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
        operations = cigar_pattern.findall(self.cigar)
        
        if not operations:
            return False
        
        # Check 5' end (first operation)
        first_len, first_op = int(operations[0][0]), operations[0][1]
        if first_op in ['S', 'H'] and first_len > max_clip_size:
            return True
        
        # Check 3' end (last operation)
        last_len, last_op = int(operations[-1][0]), operations[-1][1]
        if last_op in ['S', 'H'] and last_len > max_clip_size:
            return True
        
        return False
    
    def get_splice_junctions(self) -> Optional[Tuple[int, ...]]:
        """
        Extract splice junctions from CIGAR string
        Returns tuple of junction coordinates or None for single-block reads
        """
        # Parse CIGAR to get block structure
        cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
        operations = cigar_pattern.findall(self.cigar)
        
        if not operations:
            return None
        
        # Track genomic position and find splice junctions (N operations)
        junctions = []
        current_pos = self.start  # 0-based
        has_intron = False
        
        for length, op in operations:
            length = int(length)
            
            if op == 'M' or op == '=' or op == 'X':  # Match/mismatch
                current_pos += length
            elif op == 'N':  # Intron (splice junction)
                # Junction coordinates: last base of previous exon, first base of next exon
                # current_pos is 0-based position after last exon base
                # In 1-based: last exon base = current_pos (0-based end = 1-based last position)
                junctions.append(current_pos)      # Donor: last base of exon (1-based)
                current_pos += length
                # current_pos is now 0-based position of first base of next exon
                # In 1-based: first exon base = current_pos + 1
                junctions.append(current_pos + 1)  # Acceptor: first base of exon (1-based)
                has_intron = True
            elif op == 'D':  # Deletion
                current_pos += length
            elif op in ['I', 'S', 'H', 'P']:  # Insertion, soft/hard clip, padding
                pass  # Don't advance genomic position
        
        if not has_intron:
            return None
        
        return tuple(junctions)
    

class AlignmentProcessor:
    """Process alignments and assign reads to transcripts"""
    
    def __init__(self,
                 multi_exon_dict: Dict,
                 single_exon_dict: Dict,
                 min_mapq: int = 1,
                 min_overlap_ratio: float = 0.2):
        """
        Initialize alignment processor
        
        Args:
            multi_exon_dict: Multi-exon transcript dictionary from GTFProcessor
            single_exon_dict: Single exon dictionary from GTFProcessor
            min_mapq: Minimum mapping quality
            min_overlap_ratio: Minimum overlap ratio (intersection/transcript_length) for single-exon assignment
        """
        self.multi_exon_dict = multi_exon_dict
        self.single_exon_dict = single_exon_dict
        self.min_mapq = min_mapq
        self.min_overlap_ratio = min_overlap_ratio
        self.stats = defaultdict(int)
    
    def process_alignment(self,
                         genome_fasta: str,
                         fastq_files: List[str],
                         bed_file: str,
                         threads: int = 8) -> Dict[str, List[int]]:
        """
        Run minimap2 and process alignments on-the-fly.
        Automatically splits input files into PacBio (filename contains 'pacbio')
        and Nanopore groups, running a separate minimap2 invocation for each.
        
        Args:
            genome_fasta: Path to genome FASTA
            fastq_files: List of FASTQ files
            bed_file: Path to reference BED file (converted from GTF)
            threads: Number of threads
            
        Returns:
            Dictionary mapping transcript_id to list of read_end_positions
        """
        # Split input files by platform
        pacbio_files = [f for f in fastq_files if 'pacbio' in f.lower() or 'subread' in f.lower()]
        nanopore_files = [f for f in fastq_files if 'pacbio' not in f.lower() and 'subread' not in f.lower()]
        
        logger.info(f"Input files — Nanopore: {len(nanopore_files)}, PacBio: {len(pacbio_files)}")
        
        transcript_reads: Dict[str, list] = defaultdict(list)
        
        if nanopore_files:
            logger.info(f"Running minimap2 for Nanopore reads ({len(nanopore_files)} file(s))")
            self._run_minimap2(
                genome_fasta=genome_fasta,
                fastq_files=nanopore_files,
                bed_file=bed_file,
                threads=threads,
                preset='splice',
                extra_args=['-u', 'f', '-k', '14', '-G', '500000'],
                desc="Processing Nanopore reads",
                transcript_reads=transcript_reads,
            )
        
        if pacbio_files:
            logger.info(f"Running minimap2 for PacBio reads ({len(pacbio_files)} file(s))")
            self._run_minimap2(
                genome_fasta=genome_fasta,
                fastq_files=pacbio_files,
                bed_file=bed_file,
                threads=threads,
                preset='splice:hq',
                extra_args=['-u', 'f', '-G', '500000'],
                desc="Processing PacBio reads",
                transcript_reads=transcript_reads,
            )
        
        self._log_stats()
        logger.info(f"Assigned reads to {len(transcript_reads)} transcripts")
        
        # Sort read positions within each transcript for reproducibility
        # This ensures deterministic downstream processing regardless of alignment order
        sorted_transcript_reads = {
            transcript_id: sorted(positions)
            for transcript_id, positions in transcript_reads.items()
        }
        
        return sorted_transcript_reads

    def _run_minimap2(self,
                      genome_fasta: str,
                      fastq_files: List[str],
                      bed_file: str,
                      threads: int,
                      preset: str,
                      extra_args: List[str],
                      desc: str,
                      transcript_reads: Dict) -> None:
        """
        Run a single minimap2 invocation and process SAM output on-the-fly,
        accumulating results into the shared transcript_reads dict.
        
        Args:
            genome_fasta: Path to genome FASTA
            fastq_files: FASTQ files to align
            bed_file: BED junction guide file (may be empty string / None)
            threads: Number of threads for minimap2
            preset: minimap2 preset string (e.g. 'splice' or 'splice:hq')
            extra_args: Additional minimap2 flags inserted before reference/query
            desc: Label for the tqdm progress bar
            transcript_reads: Shared dict that results are written into
        """
        cmd = ['minimap2', '-ax', preset, '-t', str(threads)] + extra_args + ['--secondary=no']
        
        if bed_file:
            cmd.extend(['--junc-bed', bed_file])
            logger.info(f"Using junction guide from: {bed_file}")
        
        cmd.append(genome_fasta)
        cmd.extend(fastq_files)
        
        logger.info(f"Command: {' '.join(cmd)}")
        
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False,  # Binary mode for pysam
        )
        
        with pysam.AlignmentFile(process.stdout, 'r') as sam:
            for read in tqdm(sam, desc=desc):
                # Filter by mapping quality
                if read.mapping_quality < self.min_mapq:
                    self.stats['low_mapq'] += 1
                    continue
                
                # Skip unmapped, secondary, supplementary
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    self.stats['filtered_alignments'] += 1
                    continue
                
                self.stats['total_reads_processed'] += 1
                
                aligned_read = AlignedRead(
                    read_id=read.query_name,
                    chromosome=read.reference_name,
                    start=read.reference_start,  # 0-based
                    end=read.reference_end,       # 0-based
                    strand='-' if read.is_reverse else '+',
                    cigar=read.cigarstring,
                    mapq=read.mapping_quality,
                )
                
                # Filter reads with large soft/hard clips
                # if aligned_read.has_large_clips():
                #     self.stats['large_clips_filtered'] += 1
                #     continue
                
                assigned_transcript = self._assign_read(aligned_read)
                
                if not assigned_transcript:
                    self.stats['unassigned_reads'] += 1
                    continue
                
                # Determine 3' end position (1-based) based on strand
                if aligned_read.strand == '+':
                    read_end = aligned_read.end
                else:
                    read_end = aligned_read.start + 1
                
                transcript_reads[assigned_transcript].append(read_end)
                self.stats['assigned_reads'] += 1
        
        _, stderr = process.communicate()
        if process.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {stderr.decode()}")

    def _assign_read(self, read: AlignedRead) -> Optional[str]:
        """
        Assign a single read to transcript(s)
        
        Returns:
            List of assigned transcript IDs
        """
        chrom_strand_key = (read.chromosome, read.strand)
        
        # Try multi-exon assignment first
        junctions = read.get_splice_junctions()
        
        if junctions:
            # Multi-exon read - match based on splice pattern
            if chrom_strand_key in self.multi_exon_dict:
                chrom_strand_multi_exon_dict = self.multi_exon_dict[chrom_strand_key]
                
                if junctions in chrom_strand_multi_exon_dict:
                    # Single value: (transcript_id, tx_start, tx_end)
                    transcript_id, tx_start, tx_end = chrom_strand_multi_exon_dict[junctions]
                    self.stats['multi_exon_assigned'] += 1
                    return transcript_id
                else:
                    self.stats['multi_exon_no_match'] += 1
            else:
                self.stats['multi_exon_no_match'] += 1
        else:
            # Single-exon read - use InterLap for overlap-based assignment
            if chrom_strand_key in self.single_exon_dict:
                interlap = self.single_exon_dict[chrom_strand_key]
                
                # Query interval (convert to 1-based)
                read_start = read.start + 1
                read_end = read.end
                
                # Find overlapping intervals
                overlaps = list(interlap.find((read_start, read_end)))
                
                if overlaps:
                    # Find transcript with best overlap ratio
                    # overlaps contain tuples: (tx_start, tx_end, transcript_id)
                    best_transcript = None
                    best_ratio = 0.0
                    
                    for tx_start, tx_end, transcript_id in overlaps:
                        # Calculate intersection
                        intersect_start = max(read_start, tx_start)
                        intersect_end = min(read_end, tx_end)
                        intersect_length = max(0, intersect_end - intersect_start + 1)
                        
                        # Calculate transcript length
                        transcript_length = tx_end - tx_start + 1
                        
                        # Calculate overlap ratio
                        overlap_ratio = intersect_length / transcript_length
                        
                        # Track best overlap
                        if overlap_ratio > best_ratio:
                            best_ratio = overlap_ratio
                            best_transcript = transcript_id
                    
                    # Return best transcript if it meets minimum threshold
                    if best_ratio >= self.min_overlap_ratio:
                        self.stats['single_exon_assigned'] += 1
                        return best_transcript
                
                self.stats['single_exon_no_overlap'] += 1
        
        return None
    
    def _log_stats(self):
        """Log processing statistics"""
        logger.info("Read Assignment Statistics:")
        logger.info(f"  Total reads processed: {self.stats['total_reads_processed']}")
        logger.info(f"  Assigned reads: {self.stats['assigned_reads']}")
        logger.info(f"    Multi-exon: {self.stats['multi_exon_assigned']}")
        logger.info(f"    Single-exon: {self.stats['single_exon_assigned']}")
        logger.info(f"  Unassigned reads: {self.stats['unassigned_reads']}")
        logger.info(f"  Filtered reads:")
        logger.info(f"    Low MAPQ: {self.stats['low_mapq']}")
        logger.info(f"    Large clips (>100bp): {self.stats['large_clips_filtered']}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
