#!/usr/bin/env python3
"""
TES Finder: Identify alternative transcription end sites (polyA sites)
using Gaussian Mixture Models
"""

import numpy as np
import warnings
import logging
import os
import random
from typing import List, Tuple, Dict
from collections import Counter
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from dataclasses import dataclass
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Set random seeds for reproducibility
np.random.seed(42)
random.seed(42)
# Ensure single-threaded operations for full reproducibility
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'


@dataclass
class APASite:
    """Represents an alternative polyadenylation site"""
    position: int       # 1-based genomic position
    read_count: int     # Number of supporting reads


@dataclass
class TranscriptAPA:
    """APA information for a transcript"""
    site: List[int]          # List of positions (1-based)
    count: List[int]         # List of read counts per site
    abundance: List[float]   # List of relative abundances (count/total)


class TESAnalyzer:
    """Analyze termination sites using Gaussian Mixture Models"""
    
    def __init__(self,
                 min_cluster_size: int = 10,
                 max_k: int = 5,
                 min_distance: int = 50,
                 min_relative_dominance: float = 0.1,
                 min_sharpness: float = 0.5,
                 random_seed: int = 42):
        """
        Initialize TES analyzer
        
        Args:
            min_cluster_size: Minimum reads in a TES cluster
            max_k: Maximum number of GMM components to test
            min_distance: Minimum distance between TES peaks (bp)
            min_relative_dominance: Minimum relative size for alternative TES
            min_sharpness: Minimum sharpness score for peaks
            random_seed: Random seed for reproducibility
        """
        self.min_cluster_size = min_cluster_size
        self.max_k = max_k
        self.min_distance = min_distance
        self.min_relative_dominance = min_relative_dominance
        self.min_sharpness = min_sharpness
        self.random_seed = random_seed
    
    def find_tes_peaks(self, read_end_sites: np.ndarray) -> List[Tuple[int, int]]:
        """
        Find TES peaks from read end sites using all reads
        
        Args:
            read_end_sites: Array of 1-based read end positions
            
        Returns:
            List of (position, read_count) tuples for significant peaks
        """
        if len(read_end_sites) < self.min_cluster_size:
            return []
        
        # Sort read end sites for reproducibility
        # This ensures GMM receives data in consistent order regardless of
        # how reads were collected or processed
        read_end_sites = np.sort(read_end_sites)
        
        return self._find_peaks_single(read_end_sites)
    
    def _find_peaks_single(self, read_end_sites: np.ndarray) -> List[Tuple[int, int]]:
        """
        Single-pass GMM peak detection
        """
        # Find optimal number of clusters
        k_optimal = self._find_optimal_k(read_end_sites)
        
        if k_optimal == 1:
            # Single cluster - return mode
            # Handle ties deterministically by choosing smallest position
            position_counts = Counter(read_end_sites)
            max_count = max(position_counts.values())
            mode_pos = min([pos for pos, cnt in position_counts.items() if cnt == max_count])
            return [(mode_pos, max_count)]
        
        # Fit GMM with optimal k (GMM needs 2D array)
        X = read_end_sites.reshape(-1, 1)
        gmm = GaussianMixture(n_components=k_optimal, random_state=self.random_seed,
                              max_iter=200, n_init=3)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            gmm.fit(X)
        labels = gmm.predict(X)
        
        # Extract significant clusters
        return self._extract_peaks(read_end_sites, labels)
    
    def _find_optimal_k(self, read_end_sites: np.ndarray) -> int:
        """Find optimal number of clusters using silhouette score"""
        max_sil = -1
        k_optimal = 1
        n_samples = len(read_end_sites)
        X = read_end_sites.reshape(-1, 1)
        
        with warnings.catch_warnings(record=True):
            warnings.filterwarnings("ignore")
            
            for k in range(2, min(self.max_k + 1, n_samples)):
                try:
                    gmm = GaussianMixture(n_components=k, random_state=self.random_seed,
                                         max_iter=200, n_init=3)
                    gmm.fit(X)
                    labels = gmm.predict(X)
                    
                    if len(np.unique(labels)) < 2:
                        continue
                    
                    sil_score = silhouette_score(X, labels, metric='euclidean')
                    
                    if sil_score > max_sil:
                        max_sil = sil_score
                        k_optimal = k
                
                except Exception:
                    break
        
        return k_optimal
    
    def _is_peak_sharp(self, position: int, read_end_sites: np.ndarray) -> bool:
        """Check if peak is sharp using IQR method"""
        # Get reads near this peak
        window_reads = read_end_sites[
            (read_end_sites >= position - self.min_distance / 2) &
            (read_end_sites <= position + self.min_distance / 2)
        ]
        
        if len(window_reads) < self.min_cluster_size:
            return False
        
        # Calculate sharpness using IQR
        q75, q25 = np.percentile(window_reads, [75, 25])
        iqr = q75 - q25
        
        sharpness = 1.0 - (iqr / self.min_distance)
        return bool(sharpness >= self.min_sharpness)
    
    def _extract_peaks(self, read_end_sites: np.ndarray, labels: np.ndarray) -> List[Tuple[int, int]]:
        """Extract significant peaks from clustered data"""
        peaks = []
        
        # Count reads per cluster
        cluster_counts = Counter(labels)
        
        # Filter by minimum size and sort by size (descending), then by cluster_id for determinism
        # This ensures consistent ordering when clusters have equal sizes
        valid_clusters = sorted(
            [(k, v) for k, v in cluster_counts.items() if v >= self.min_cluster_size],
            key=lambda x: (x[1], x[0]), reverse=True
        )
        
        if not valid_clusters:
            return []
        
        dominant_size = valid_clusters[0][1]
        used_positions = set()
        
        for cluster_id, cluster_size in valid_clusters:
            # Get positions for this cluster
            cluster_positions = read_end_sites[labels == cluster_id]
            
            # Find mode position (most common position)
            # If there are ties, Counter.most_common() returns them in arbitrary order,
            # so we need to handle ties deterministically by choosing the smallest position
            position_counts = Counter(cluster_positions)
            max_count = max(position_counts.values())
            # Get all positions with max count and choose the smallest for determinism
            mode_position = min([pos for pos, cnt in position_counts.items() if cnt == max_count])
            
            # Check sharpness
            if not self._is_peak_sharp(mode_position, read_end_sites):
                continue
            
            # Check relative dominance
            if cluster_size < self.min_relative_dominance * dominant_size:
                continue
            
            # Check distance from other peaks
            too_close = any(abs(mode_position - pos) < self.min_distance
                           for pos in used_positions)
            
            if not too_close:
                peaks.append((mode_position, cluster_size))
                used_positions.add(mode_position)
        
        return sorted(peaks, key=lambda x: x[0])


class TESFinder:
    """Parallel TES Finder for identifying alternative polyA sites"""
    
    def __init__(self,
                 min_reads: int = 10,
                 min_cluster_size: int = 10,
                 max_k: int = 5,
                 min_distance: int = 50,
                 min_relative_dominance: float = 0.1,
                 min_sharpness: float = 0.5,
                 n_jobs: int = -1,
                 random_seed: int = 42):
        """
        Initialize TES Finder
        
        Args:
            min_reads: Minimum reads to analyze a transcript
            min_cluster_size: Minimum reads in a TES cluster
            max_k: Maximum GMM components
            min_distance: Minimum distance between peaks (bp)
            min_relative_dominance: Minimum relative size for alternative site
            min_sharpness: Minimum peak sharpness
            n_jobs: Number of parallel jobs (-1 for all CPUs, 1 for serial/reproducible)
            random_seed: Random seed for reproducibility
        """
        self.min_reads = min_reads
        self.n_jobs = n_jobs
        self.random_seed = random_seed
        
        self.analyzer_params = {
            'min_cluster_size': min_cluster_size,
            'max_k': max_k,
            'min_distance': min_distance,
            'min_relative_dominance': min_relative_dominance,
            'min_sharpness': min_sharpness,
            'random_seed': random_seed
        }
        
        logger.info(f"Initialized TES Finder with {n_jobs} workers (random_seed={random_seed})")
    
    def find_apa_sites(self,
                       transcript_reads: Dict[str, List[int]]) -> Dict[str, TranscriptAPA]:
        """
        Find alternative polyA sites for all transcripts
        
        Args:
            transcript_reads: Dict mapping transcript_id to list of read_end_positions
            
        Returns:
            Dictionary of transcript_id -> TranscriptAPA
        """
        logger.info(f"Finding APA sites for {len(transcript_reads)} transcripts")
        
        # Prepare data for parallel processing
        # Sort positions within each transcript for reproducibility
        valid_transcripts = [
            (tid, np.sort(np.array(positions)))
            for tid, positions in transcript_reads.items()
            if len(positions) >= self.min_reads
        ]
        
        logger.info(f"Analyzing {len(valid_transcripts)} transcripts with sufficient reads")
        
        # Choose backend for reproducibility
        # Note: 'loky' is faster but may introduce non-determinism
        # Use 'threading' or set n_jobs=1 for full reproducibility
        backend = 'loky' if self.n_jobs != 1 else 'sequential'
        if self.n_jobs == 1:
            logger.info("Running in sequential mode for full reproducibility")
        
        # Process in parallel (or sequentially if n_jobs=1)
        results = Parallel(n_jobs=self.n_jobs, backend=backend)(
            delayed(self._process_transcript)(
                transcript_id,
                end_positions,
                self.analyzer_params
            ) for transcript_id, end_positions in tqdm(valid_transcripts, desc="TES Analysis")
        )

        # Shut down loky worker pool immediately so the process does not hang
        # waiting for the default pool reuse timeout (~300 s)
        if self.n_jobs != 1:
            get_reusable_executor().shutdown(wait=True)

        # Collect results
        apa_results = {}
        n_apa = 0
        
        for transcript_id, apa_info in results:
            apa_results[transcript_id] = apa_info
            if len(apa_info.site) > 1:
                n_apa += 1
        
        logger.info(f"Found alternative TES in {n_apa} transcripts")
        
        return apa_results
    
    @staticmethod
    def _process_transcript(transcript_id: str,
                           end_positions: np.ndarray,
                           analyzer_params: Dict) -> Tuple[str, TranscriptAPA]:
        """Process a single transcript (for parallel execution)"""
        
        analyzer = TESAnalyzer(**analyzer_params)
        
        # Find TES peaks
        peaks = analyzer.find_tes_peaks(end_positions)
        
        # If no peaks found, use mode
        if not peaks:
            # Handle ties deterministically by choosing smallest position
            position_counts = Counter(end_positions)
            max_count = max(position_counts.values())
            mode_pos = min([pos for pos, cnt in position_counts.items() if cnt == max_count])
            
            return transcript_id, TranscriptAPA(
                site=[int(mode_pos)],
                count=[max_count],
                abundance=[1.0]
            )
        
        # Sort peaks by read count (descending), then by position (ascending) for determinism
        peaks = sorted(peaks, key=lambda x: (-x[1], x[0]))
        
        # Split into separate lists and calculate abundance
        sites = [pos for pos, _ in peaks]
        counts = [count for _, count in peaks]
        total_count = sum(counts)
        abundances = [count / total_count for count in counts]
        
        return transcript_id, TranscriptAPA(
            site=sites,
            count=counts,
            abundance=abundances
        )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
