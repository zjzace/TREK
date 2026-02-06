#!/usr/bin/env python3
"""
ApaFinder: Main pipeline for identifying alternative polyA sites from long-read data

This pipeline processes long-read sequencing data to identify alternative 
transcription end sites (TES) / polyadenylation sites:
1. Process GTF annotation to extract transcript structures
2. Align reads using minimap2 with splice-aware mapping
3. Assign reads to transcripts based on splice junction patterns
4. Identify alternative TES using Gaussian Mixture Models
"""

import argparse
import logging
import sys
import subprocess
import pickle
from pathlib import Path

from gtf_processor import GTFProcessor
from alignment_processor import AlignmentProcessor
from apa_finder import TESFinder
from internal_priming_filter import InternalPrimingFilter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ApaFinderPipeline:
    """Main pipeline for polyA site identification"""
    
    def __init__(self, gtf_file, genome_fasta, fastq_files, output_dir, prefix='polyA',
                 threads=8, min_mapq=1, min_overlap_ratio=0.2,
                 min_reads=10, min_cluster_size=10, max_clusters=5, min_distance=50,
                 min_relative_dominance=0.1, min_sharpness=0.5, n_jobs=-1,
                priming_window=10, priming_a_threshold=0.5,
                 random_seed=42):
        """Initialize pipeline with parameters"""
        self.gtf_file = gtf_file
        self.genome_fasta = genome_fasta
        self.fastq_files = fastq_files
        self.output_dir = Path(output_dir)
        self.prefix = prefix
        self.threads = threads
        self.min_mapq = min_mapq
        self.min_overlap_ratio = min_overlap_ratio
        self.min_reads = min_reads
        self.min_cluster_size = min_cluster_size
        self.max_clusters = max_clusters
        self.min_distance = min_distance
        self.min_relative_dominance = min_relative_dominance
        self.min_sharpness = min_sharpness
        self.n_jobs = n_jobs
        self.priming_window = priming_window
        self.priming_a_threshold = priming_a_threshold
        self.random_seed = random_seed
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("ApaFinder Pipeline Started")
        logger.info(f"Output directory: {self.output_dir}")
    
    def run(self):
        """Run the complete pipeline"""
        try:
            logger.info("STEP 1: Processing GTF annotation")
            transcripts, multi_exon_dict, single_exon_dict = self._process_gtf()
            
            logger.info("STEP 2: Converting GTF to BED for alignment guide")
            bed_file = self._gtf_to_bed()
            
            logger.info("STEP 3: Running alignment and assigning reads to transcripts")
            transcript_reads = self._process_alignment(bed_file, multi_exon_dict, single_exon_dict)
            
            logger.info("STEP 4: Identifying alternative polyA sites")
            apa_results = self._find_apa_sites(transcript_reads)
            
            logger.info("STEP 5: Filtering internal priming artifacts")
            apa_results = self._filter_internal_priming(apa_results, transcripts)
            
            logger.info("STEP 6: Writing results")
            self._write_results(apa_results, transcripts)
            
            logger.info("Pipeline completed successfully!")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}", exc_info=True)
            raise
    
    def _process_gtf(self):
        """Process GTF file"""
        logger.info(f"Processing GTF: {self.gtf_file}")
        
        processor = GTFProcessor()
        transcripts = processor.parse_gtf(self.gtf_file)
        
        logger.info(f"Parsed {len(transcripts)} transcripts")
        
        multi_exon_dict, single_exon_dict = processor.organize_transcripts(transcripts)
        return transcripts, multi_exon_dict, single_exon_dict
    
    def _gtf_to_bed(self):
        """Convert GTF to BED for minimap2 junction guide"""
        bed_file = self.output_dir / f"{self.prefix}.junctions.bed"
        gtf2bed_script = Path(__file__).parent / 'gtf2bed.pl'
        
        if not gtf2bed_script.exists():
            logger.warning("gtf2bed.pl not found - alignment will proceed without junction guide")
            return None
        
        logger.info("Converting GTF to BED")
        with open(bed_file, 'w') as f:
            subprocess.run(['perl', str(gtf2bed_script), self.gtf_file], stdout=f, check=True)
        
        logger.info(f"Created BED file: {bed_file}")
        return str(bed_file)
    
    def _process_alignment(self, bed_file, multi_exon_dict, single_exon_dict):
        """Run alignment and assign reads on-the-fly"""
        processor = AlignmentProcessor(
            multi_exon_dict=multi_exon_dict,
            single_exon_dict=single_exon_dict,
            min_mapq=self.min_mapq,
            min_overlap_ratio=self.min_overlap_ratio
        )
        
        transcript_reads = processor.process_alignment(
            genome_fasta=self.genome_fasta,
            fastq_files=self.fastq_files,
            bed_file=bed_file,
            threads=self.threads
        )
        
        # Save read assignments as pickle
        assignment_file = self.output_dir / f"{self.prefix}.read_assignments.pkl"
        with open(assignment_file, 'wb') as f:
            pickle.dump(transcript_reads, f)
        
        logger.info(f"Saved read assignments: {assignment_file}")
        return transcript_reads
    
    def _find_apa_sites(self, transcript_reads):
        """Find alternative polyA sites"""
        finder = TESFinder(
            min_reads=self.min_reads,
            min_cluster_size=self.min_cluster_size,
            max_k=self.max_clusters,
            min_distance=self.min_distance,
            min_relative_dominance=self.min_relative_dominance,
            min_sharpness=self.min_sharpness,
            n_jobs=self.n_jobs,
            random_seed=self.random_seed
        )
        
        return finder.find_apa_sites(transcript_reads)
    
    def _filter_internal_priming(self, apa_results, transcripts):
        """Filter APA sites to remove internal priming artifacts"""
        filter_obj = InternalPrimingFilter(
            genome_fasta=self.genome_fasta,
            window_size=self.priming_window,
            a_content_threshold=self.priming_a_threshold
        )
        
        return filter_obj.filter_apa_results(apa_results, transcripts)
    
    def _write_results(self, apa_results, transcripts):
        """Write results to output files"""
        # Main results file
        results_file = self.output_dir / f"{self.prefix}.apa_sites.txt"
        
        with open(results_file, 'w') as f:
            f.write("transcript_id\tgene_id\tgene_name\tchromosome\tstrand\t"
                   "ID\tsite_position\tsite_count\tsite_abundance\n")
            
            for transcript_id, apa in apa_results.items():
                transcript = transcripts.get(transcript_id)
                if not transcript:
                    continue
                
                # Write one line per APA site
                for position, count, abundance in zip(apa.site, apa.count, apa.abundance):
                    # Create locus ID: chrom:position:strand
                    locus_id = f"{transcript.chromosome}:{position}:{transcript.strand}"
                    
                    f.write(f"{transcript_id}\t{transcript.gene_id}\t{transcript.gene_name}\t"
                           f"{transcript.chromosome}\t{transcript.strand}\t{locus_id}\t"
                           f"{position}\t{count}\t{abundance:.4f}\n")
        
        logger.info(f"Saved results: {results_file}")
        
        # Summary file
        summary_file = self.output_dir / f"{self.prefix}.summary.txt"
        
        with open(summary_file, 'w') as f:
            total = len(apa_results)
            with_apa = sum(1 for apa in apa_results.values() if len(apa.site) > 1)
            
            f.write("ApaFinder Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Total transcripts analyzed: {total}\n")
            f.write(f"Transcripts with alternative TES: {with_apa}\n")
            f.write(f"Percentage with APA: {100 * with_apa / total:.2f}%\n")
            
            site_counts = {}
            for apa in apa_results.values():
                n = len(apa.site)
                site_counts[n] = site_counts.get(n, 0) + 1
            
            f.write("\nDistribution of TES per transcript:\n")
            for n in sorted(site_counts.keys()):
                f.write(f"  {n} sites: {site_counts[n]} transcripts\n")
        
        logger.info(f"Saved summary: {summary_file}")


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='ApaFinder: Identify alternative polyA sites from long-read data'
    )
    
    # Required arguments
    parser.add_argument('-g', '--gtf', required=True, 
                       help='Reference GTF file')
    parser.add_argument('-f', '--fasta', required=True, 
                       help='Genome FASTA file')
    parser.add_argument('-q', '--fastq', nargs='+', required=True, 
                       help='Input FASTQ.gz file(s)')
    
    # Output arguments
    parser.add_argument('-o', '--output', default='results', 
                       help='Output directory (default: results)')
    parser.add_argument('-p', '--prefix', default='polyA', 
                       help='Output file prefix (default: polyA)')
    
    # Alignment arguments
    parser.add_argument('-t', '--threads', type=int, default=8, 
                       help='Number of threads (default: 8)')
    parser.add_argument('--min-mapq', type=int, default=1, 
                       help='Minimum MAPQ (default: 1)')
    parser.add_argument('--min-overlap', type=float, default=0.2, 
                       help='Minimum overlap ratio (default: 0.2)')
    
    # TES detection arguments
    parser.add_argument('--min-reads', type=int, default=10, 
                       help='Minimum reads per transcript (default: 10)')
    parser.add_argument('--min-cluster-size', type=int, default=10, 
                       help='Minimum cluster size (default: 10)')
    parser.add_argument('--max-clusters', type=int, default=5, 
                       help='Maximum clusters (default: 5)')
    parser.add_argument('--min-distance', type=int, default=50, 
                       help='Minimum distance between sites in bp (default: 50)')
    parser.add_argument('--min-dominance', type=float, default=0.3, 
                       help='Minimum relative dominance (default: 0.3)')
    parser.add_argument('--min-sharpness', type=float, default=0.5, 
                       help='Minimum peak sharpness (default: 0.5)')
    
    # Performance arguments
    parser.add_argument('-j', '--jobs', type=int, default=-1, 
                       help='Parallel jobs (default: -1, all CPUs)')
    
    # Internal priming filter arguments
    parser.add_argument('--priming-window', type=int, default=10,
                       help='Upstream window size for A-content check (default: 10)')
    parser.add_argument('--priming-a-threshold', type=float, default=0.5,
                       help='Maximum A-content threshold (default: 0.5)')
    
    # Reproducibility arguments
    parser.add_argument('--serial', action='store_true',
                       help='Run in serial mode for full reproducibility (overrides -j)')
    
    return parser.parse_args()


def main():
    """Main entry point"""
    args = parse_arguments()
    
    try:
        # Validate inputs
        if not Path(args.gtf).exists():
            raise FileNotFoundError(f"GTF file not found: {args.gtf}")
        if not Path(args.fasta).exists():
            raise FileNotFoundError(f"FASTA file not found: {args.fasta}")
        for fq in args.fastq:
            if not Path(fq).exists():
                raise FileNotFoundError(f"FASTQ file not found: {fq}")
        
        # Run pipeline
        pipeline = ApaFinderPipeline(
            gtf_file=args.gtf,
            genome_fasta=args.fasta,
            fastq_files=args.fastq,
            output_dir=args.output,
            prefix=args.prefix,
            threads=args.threads,
            min_mapq=args.min_mapq,
            min_overlap_ratio=args.min_overlap,
            min_reads=args.min_reads,
            min_cluster_size=args.min_cluster_size,
            max_clusters=args.max_clusters,
            min_distance=args.min_distance,
            min_relative_dominance=args.min_dominance,
            min_sharpness=args.min_sharpness,
            n_jobs=1 if args.serial else args.jobs,
            priming_window=args.priming_window,
            priming_a_threshold=args.priming_a_threshold,
            random_seed=42
        )
        pipeline.run()
        
        return 0
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())
