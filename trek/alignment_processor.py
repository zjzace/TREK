#!/usr/bin/env python3
"""
Alignment Processor: Handle minimap2 alignment and read-to-transcript assignment
"""

import subprocess
import logging
import pysam
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from tqdm import tqdm

logger = logging.getLogger(__name__)


# pysam CIGAR operation codes
_OP_M = 0   # Match/mismatch
_OP_D = 2   # Deletion
_OP_N = 3   # Intron (splice junction)
_OP_EQ = 7  # Sequence match
_OP_X = 8   # Sequence mismatch


def _extract_junctions(cigartuples, ref_start: int) -> Optional[Tuple[int, ...]]:
    """
    Extract splice junctions from pysam cigartuples.
    
    Args:
        cigartuples: List of (op, length) from pysam
        ref_start: 0-based reference start
        
    Returns:
        Tuple of junction coordinates (1-based) or None for single-block reads
    """
    junctions = []
    pos = ref_start  # 0-based
    for op, length in cigartuples:
        if op == _OP_M or op == _OP_EQ or op == _OP_X:
            pos += length
        elif op == _OP_N:
            junctions.append(pos)          # Donor: 0-based end = 1-based last position
            pos += length
            junctions.append(pos + 1)      # Acceptor: 1-based first position
        elif op == _OP_D:
            pos += length
        # I(1), S(4), H(5), P(6) don't advance reference position
    return tuple(junctions) if junctions else None


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
                extra_args=['-u', 'b', '-k', '14', '-G', '500000'],
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
        base_cmd = ['minimap2', '-ax', preset, '-t', str(threads)] + extra_args + ['--secondary=no']
        
        if bed_file:
            base_cmd.extend(['--junc-bed', bed_file])
            logger.info(f"Using junction guide from: {bed_file}")
        
        base_cmd.append(genome_fasta)
        
        # Process each FASTQ file independently so one truncated file
        # doesn't abort the entire batch
        for file_idx, fastq_file in enumerate(fastq_files, 1):
            cmd = base_cmd + [fastq_file]
            file_desc = f"{desc} [{file_idx}/{len(fastq_files)}] {fastq_file}"
            logger.info(f"Command: {' '.join(cmd)}")
            
            try:
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=False,  # Binary mode for pysam
                )
                
                with pysam.AlignmentFile(process.stdout, 'r') as sam:
                  try:
                    for read in tqdm(sam, desc=file_desc):
                        # Filter by mapping quality
                        if read.mapping_quality < self.min_mapq:
                            self.stats['low_mapq'] += 1
                            continue
                        
                        # Skip unmapped, secondary, supplementary
                        if read.is_unmapped or read.is_secondary or read.is_supplementary:
                            self.stats['filtered_alignments'] += 1
                            continue
                        
                        self.stats['total_reads_processed'] += 1
                        
                        chrom = read.reference_name
                        ref_start = read.reference_start   # 0-based
                        ref_end = read.reference_end         # 0-based
                        strand = '-' if read.is_reverse else '+'
                        
                        junctions = _extract_junctions(read.cigartuples, ref_start)
                        assigned = self._assign_read(chrom, strand, ref_start, ref_end, junctions)
                        
                        if not assigned:
                            self.stats['unassigned_reads'] += 1
                            continue
                        
                        assigned_transcript, transcript_strand = assigned
                        
                        # Determine 3' end position (1-based) based on transcript strand
                        if transcript_strand == '+':
                            read_end = ref_end
                        else:
                            read_end = ref_start + 1
                        
                        transcript_reads[assigned_transcript].append(read_end)
                        self.stats['assigned_reads'] += 1
                  except (StopIteration, IOError, OSError) as e:
                    logger.warning(f"SAM stream interrupted for {fastq_file} "
                                   f"(truncated file?): {e}. "
                                   f"Reads processed so far from this file are kept.")
                
                _, stderr = process.communicate()
                if process.returncode != 0:
                    logger.warning(f"minimap2 returned non-zero exit for {fastq_file}: "
                                   f"{stderr.decode().strip()}. Continuing with next file.")
                    self.stats['failed_files'] += 1
            
            except Exception as e:
                logger.warning(f"Failed to process {fastq_file}: {e}. Skipping.")
                self.stats['failed_files'] += 1

    def _assign_read(self, chrom: str, strand: str,
                     ref_start: int, ref_end: int,
                     junctions: Optional[Tuple[int, ...]]) -> Optional[Tuple[str, str]]:
        """
        Assign a single read to a transcript
        
        Args:
            chrom: Chromosome name
            strand: Read strand ('+' or '-')
            ref_start: 0-based reference start
            ref_end: 0-based reference end
            junctions: Splice junction tuple from _extract_junctions, or None
        
        Returns:
            Tuple of (transcript_id, transcript_strand), or None if unassigned
        """
        if junctions:
            # Multi-exon read - match based on splice pattern (strand-agnostic)
            chrom_multi_exon_dict = self.multi_exon_dict.get(chrom)
            if chrom_multi_exon_dict:
                match = chrom_multi_exon_dict.get(junctions)
                if match:
                    transcript_id, tx_start, tx_end, tx_strand = match
                    self.stats['multi_exon_assigned'] += 1
                    return transcript_id, tx_strand
                else:
                    self.stats['multi_exon_no_match'] += 1
            else:
                self.stats['multi_exon_no_match'] += 1
        else:
            # Single-exon read - use InterLap for overlap-based assignment
            chrom_strand_key = (chrom, strand)
            interlap = self.single_exon_dict.get(chrom_strand_key)
            if interlap:
                # Query interval (convert to 1-based)
                read_start = ref_start + 1
                read_end = ref_end
                
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
                        return best_transcript, strand
                
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


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
