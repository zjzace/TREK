#!/usr/bin/env python3
"""
Internal Priming Filter: Remove APA sites that may result from internal priming 
by checking A-rich sequences around the polyA site
"""

import logging
from typing import Dict
from Bio import SeqIO
from apa_finder import TranscriptAPA
from gtf_processor import Transcript

logger = logging.getLogger(__name__)


class InternalPrimingFilter:
    """Filter to remove APA sites that may result from internal priming"""
    
    def __init__(self, 
                 genome_fasta: str,
                 window_size: int = 10,
                 a_content_threshold: float = 0.5):
        """
        Initialize internal priming filter
        
        Args:
            genome_fasta: Path to genome FASTA file
            window_size: Size of upstream window to check (bp)
            a_content_threshold: Maximum allowed A content (0-1)
        """
        self.genome_fasta = genome_fasta
        self.window_size = window_size
        self.a_content_threshold = a_content_threshold
        
        # Load genome sequences
        logger.info(f"Loading genome FASTA: {genome_fasta}")
        self.genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
        logger.info(f"Loaded {len(self.genome_seqs)} chromosomes")
    
    def filter_apa_results(self,
                          apa_results: Dict[str, TranscriptAPA],
                          transcripts: Dict[str, Transcript]) -> Dict[str, TranscriptAPA]:
        """
        Filter APA results to remove sites with internal priming
        
        Args:
            apa_results: Dictionary of transcript_id -> TranscriptAPA
            transcripts: Dictionary of transcript_id -> Transcript
            
        Returns:
            Filtered dictionary of transcript_id -> TranscriptAPA
        """
        logger.info("Filtering APA sites for internal priming")
        
        filtered_results = {}
        total_sites_before = 0
        total_sites_after = 0
        transcripts_with_multiple_before = 0
        transcripts_with_multiple_after = 0
        sites_removed = 0
        
        for transcript_id, apa in apa_results.items():
            transcript = transcripts.get(transcript_id)
            if not transcript:
                logger.warning(f"Transcript {transcript_id} not found in GTF")
                filtered_results[transcript_id] = apa
                continue
            
            total_sites_before += len(apa.site)
            
            # Only filter transcripts with multiple APA sites
            if len(apa.site) <= 1:
                filtered_results[transcript_id] = apa
                total_sites_after += len(apa.site)
            else:
                transcripts_with_multiple_before += 1
                
                # Filter sites
                filtered_apa = self._filter_transcript_apa(apa, transcript)
                filtered_results[transcript_id] = filtered_apa
                
                total_sites_after += len(filtered_apa.site)
                sites_removed += (len(apa.site) - len(filtered_apa.site))
                
                if len(filtered_apa.site) > 1:
                    transcripts_with_multiple_after += 1
        
        logger.info(f"Internal priming filter results:")
        logger.info(f"  Total APA sites before: {total_sites_before}")
        logger.info(f"  Total APA sites after: {total_sites_after}")
        logger.info(f"  Sites removed: {sites_removed}")
        logger.info(f"  Transcripts with multiple sites before: {transcripts_with_multiple_before}")
        logger.info(f"  Transcripts with multiple sites after: {transcripts_with_multiple_after}")
        
        return filtered_results
    
    def _filter_transcript_apa(self,
                               apa: TranscriptAPA,
                               transcript: Transcript) -> TranscriptAPA:
        """
        Filter APA sites for a single transcript
        Only filters non-dominant sites; the dominant site (highest abundance) is always kept.
        
        Args:
            apa: TranscriptAPA object (sites sorted by count, descending)
            transcript: Transcript object
            
        Returns:
            Filtered TranscriptAPA object
        """
        # Keep track of which sites to keep
        # Always keep the dominant site (index 0)
        keep_indices = [0]
        
        # Only check non-dominant sites (index >= 1)
        for idx in range(1, len(apa.site)):
            position = apa.site[idx]
            
            # Extract sequence around this position
            a_proportion = self._calculate_a_content(
                chromosome=transcript.chromosome,
                position=position,
                strand=transcript.strand
            )
            
            # Keep site if A proportion is below threshold
            if a_proportion <= self.a_content_threshold:
                keep_indices.append(idx)
                logger.debug(f"Keeping non-dominant site at {position} in {transcript.transcript_id} "
                           f"(A proportion: {a_proportion:.2%})")
            else:
                logger.debug(f"Removing non-dominant site at {position} in {transcript.transcript_id} "
                           f"(A proportion: {a_proportion:.2%})")
        
        # Create filtered APA object
        filtered_sites = [apa.site[i] for i in keep_indices]
        filtered_counts = [apa.count[i] for i in keep_indices]
        
        # Recalculate abundances
        total_count = sum(filtered_counts)
        filtered_abundances = [count / total_count for count in filtered_counts]
        
        return TranscriptAPA(
            site=filtered_sites,
            count=filtered_counts,
            abundance=filtered_abundances
        )
    
    def _calculate_a_content(self,
                            chromosome: str,
                            position: int,
                            strand: str) -> float:
        """
        Calculate A proportion in upstream window of APA site
        
        Args:
            chromosome: Chromosome name
            position: Genomic position (1-based)
            strand: Strand (+ or -)
            
        Returns:
            A proportion as fraction (0-1)
        """
        # Get chromosome sequence
        if chromosome not in self.genome_seqs:
            logger.warning(f"Chromosome {chromosome} not found in genome")
            return 0.0
        
        chr_seq = self.genome_seqs[chromosome].seq
        
        # Calculate upstream window coordinates
        # Position is 1-based, convert to 0-based for sequence extraction
        if strand == '+':
            # For + strand: upstream means lower genomic coordinates
            # Extract window_size bp immediately upstream of the APA site
            # [position - window_size, position - 1] in 1-based
            # [position - window_size - 1, position - 1] in 0-based
            start = max(0, position - self.window_size - 1)
            end = max(0, position - 1)
        else:
            # For - strand: upstream in transcript direction means higher genomic coordinates
            # Extract window_size bp immediately after the APA site
            # [position + 1, position + window_size] in 1-based
            # [position, position + window_size] in 0-based
            start = position
            end = min(len(chr_seq), position + self.window_size)
        
        # Extract sequence
        window_seq = str(chr_seq[start:end]).upper()
        
        if not window_seq:
            return 0.0
        
        # Calculate A proportion
        # For + strand: count A directly
        # For - strand: count T (which represents A in the transcript)
        if strand == '+':
            a_count = window_seq.count('A')
        else:
            a_count = window_seq.count('T')  # T in genomic = A in transcript
        
        a_proportion = a_count / len(window_seq) if len(window_seq) > 0 else 0.0
        
        return a_proportion
