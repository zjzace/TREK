# polyAFinder

A professional tool for identifying alternative polyadenylation (APA) sites from long-read sequencing data.

## Overview

polyAFinder analyzes long-read RNA sequencing data (e.g., Oxford Nanopore, PacBio) to identify alternative transcription end sites (TES) and polyadenylation sites. The tool uses Gaussian Mixture Models (GMM) to cluster read end positions and identify statistically significant alternative polyA sites.

## Features

- **Splice-aware read alignment** using minimap2 with junction guidance
- **Accurate read-to-transcript assignment** based on splice junction patterns
- **Robust TES detection** using Gaussian Mixture Models with consensus calling
- **Parallel processing** for efficient analysis of large datasets
- **Comprehensive output** with detailed APA site information

## Pipeline Overview

```
1. GTF Processing
   ├─ Parse reference GTF
   ├─ Extract transcript structures
   └─ Build splice junction patterns

2. Read Alignment (On-the-fly)
   ├─ Convert GTF to BED for junction guide
   ├─ Run minimap2 with splice-aware mapping
   ├─ Process alignments on-the-fly (no BAM output)
   ├─ Assign reads to transcripts
   └─ Collect read 3' end positions

3. TES Detection
   ├─ Cluster read end sites using GMM
   ├─ Apply consensus calling for large datasets
   ├─ Filter peaks by sharpness and dominance
   └─ Identify alternative polyA sites

4. Output
   ├─ APA site positions and read counts
   ├─ Transcript-level statistics
   └─ Summary report
```

## Installation

### Requirements

- Python >= 3.7
- minimap2
- samtools
- Perl (for GTF to BED conversion)

### Python Dependencies

Install required Python packages:

```bash
pip install numpy pandas scikit-learn pysam PyYAML joblib tqdm
```

Or use the requirements file:

```bash
pip install -r requirements.txt
```

### External Tools

Install minimap2 and samtools:

```bash
# On Ubuntu/Debian
sudo apt-get install minimap2 samtools

# On macOS with Homebrew
brew install minimap2 samtools

# Using conda
conda install -c bioconda minimap2 samtools
```

## Quick Start

### 1. Prepare Configuration File

Create a `config.yaml` file (see `config.yaml` template):

```yaml
input:
  reference_gtf: "reference.gtf"
  genome_fasta: "genome.fa"
  fastq_files:
    - "sample1.fastq.gz"
    - "sample2.fastq.gz"

output:
  output_dir: "results"
  prefix: "my_sample"

alignment:
  threads: 8

tes_detection:
  min_reads: 10
  min_cluster_size: 10
  min_distance: 50
```

### 2. Run the Pipeline

```bash
# Using configuration file
python polyAFinder/polyAFinder.py -c config.yaml

# Using command line arguments
python polyAFinder/polyAFinder.py \
  -g reference.gtf \
  -f genome.fa \
  -q reads1.fastq.gz reads2.fastq.gz \
  -o results \
  -t 8
```

## Usage

### Command Line Options

```
polyAFinder.py [-h] [-c CONFIG] [-g GTF] [-f FASTA] [-q FASTQ [FASTQ ...]]
               [-o OUTPUT] [-p PREFIX] [-t THREADS] [-j JOBS]

Arguments:
  -c, --config CONFIG       Configuration YAML file
  -g, --gtf GTF            Reference GTF file
  -f, --fasta FASTA        Genome FASTA file
  -q, --fastq FASTQ        Input FASTQ.gz file(s)
  -o, --output OUTPUT      Output directory (default: results)
  -p, --prefix PREFIX      Output file prefix (default: polyA)
  -t, --threads THREADS    Number of alignment threads (default: 8)
  -j, --jobs JOBS          Number of parallel jobs for analysis (default: -1)
```

### Configuration Parameters

#### Input/Output
- `reference_gtf`: Path to reference GTF annotation
- `genome_fasta`: Path to genome FASTA file
- `fastq_files`: List of input FASTQ files
- `output_dir`: Output directory path
- `prefix`: Prefix for output files

#### Alignment Settings
- `threads`: Number of threads for minimap2
- `min_mapq`: Minimum mapping quality (default: 10)

#### GTF Processing
- `overlap_tolerance`: Tolerance for single-exon overlap (default: 20)

#### Read Assignment
- `perfect_match_only`: Only assign multi-exon reads with perfect SJ match (default: true)

#### TES Detection
- `min_reads`: Minimum reads supporting a transcript (default: 10)
- `min_cluster_size`: Minimum reads in a TES cluster (default: 10)
- `max_clusters`: Maximum GMM components to test (default: 5)
- `min_distance`: Minimum distance between TES peaks in bp (default: 50)
- `min_relative_dominance`: Minimum relative size for alternative TES (default: 0.1)
- `min_sharpness`: Minimum sharpness score for peaks (default: 0.5)
- `downsample_threshold`: Threshold for downsampling (default: 2000)
- `n_iterations`: Number of consensus iterations (default: 5)

#### Performance
- `n_jobs`: Number of parallel jobs (-1 = all CPUs)

## Output Files

### 1. `{prefix}.apa_sites.txt`

Main results file with APA site information:

| Column | Description |
|--------|-------------|
| transcript_id | Transcript identifier |
| gene_id | Gene identifier |
| chromosome | Chromosome name |
| strand | Strand (+/-) |
| has_apa | Whether alternative TES detected |
| num_sites | Number of TES sites |
| dominant_position | Position of dominant TES |
| site_positions | Comma-separated positions of all sites |
| site_counts | Comma-separated read counts for each site |

### 2. `{prefix}.read_assignments.txt`

Read-to-transcript assignments:

| Column | Description |
|--------|-------------|
| transcript_id | Assigned transcript |
| end_position | Read 3' end position (0-based) |

### 3. `{prefix}.summary.txt`

Summary statistics:
- Total transcripts analyzed
- Transcripts with alternative TES
- Distribution of TES per transcript

### 4. `{prefix}.log`

Detailed log file with pipeline execution information

## Algorithm Details

### TES Detection Using GMM

polyAFinder uses a sophisticated Gaussian Mixture Model approach to identify alternative polyA sites:

1. **Small datasets** (≤2000 reads): Single-pass GMM clustering
2. **Large datasets** (>2000 reads): Consensus approach with multiple iterations

For each dataset:
- Determine optimal number of clusters using silhouette score
- Fit GMM to cluster read end positions
- Filter clusters by:
  - Minimum size
  - Peak sharpness (using IQR)
  - Relative dominance
  - Minimum distance between peaks

### Read Assignment Strategy

**Multi-exon reads:**
- Extract splice junction pattern from CIGAR
- Match against reference transcript splice patterns
- Require perfect junction match
- Verify chromosome, strand, and coordinate consistency

**Single-exon reads:**
- Find overlapping single-exon transcripts
- Allow configurable overlap tolerance
- Match chromosome and strand

## Examples

### Example 1: Basic Usage

```bash
python polyAFinder/polyAFinder.py \
  -g data/reference.gtf \
  -f data/genome.fa \
  -q data/sample.fastq.gz \
  -o results/sample1 \
  -t 16
```

### Example 2: Multiple FASTQ Files

```bash
python polyAFinder/polyAFinder.py \
  -g data/reference.gtf \
  -f data/genome.fa \
  -q data/rep1.fastq.gz data/rep2.fastq.gz data/rep3.fastq.gz \
  -o results/merged \
  -t 32 \
  -j 16
```

### Example 3: Using Configuration File

```bash
python polyAFinder/polyAFinder.py -c my_config.yaml
```

## Performance Considerations

- **Memory**: Moderate RAM usage (~2-4 GB for typical datasets) - no intermediate BAM files
- **Threads**: Use `-t` for alignment and processing threads, `-j` for TES analysis parallelization
- **On-the-fly processing**: Alignments are processed directly from minimap2 output without disk I/O
- **Large datasets**: Automatic downsampling with consensus calling ensures robustness
- **Disk space**: Minimal - only output files are written (no intermediate BAM files)

## Troubleshooting

### Common Issues

1. **"GTF file not found"**: Verify file paths in configuration
2. **"minimap2: command not found"**: Install minimap2 and ensure it's in PATH
3. **Low assignment rate**: Check GTF compatibility with genome version
4. **Memory errors**: Reduce `-j` parameter or increase system memory

### Getting Help

For issues or questions:
1. Check the log file for detailed error messages
2. Verify input file formats (GTF, FASTA, FASTQ)
3. Ensure all dependencies are installed

## Citation

If you use polyAFinder in your research, please cite:

```
[Your citation information here]
```

## License

This project is licensed under the MIT License.

## Acknowledgments

This tool was developed based on methods from:
- LAFITE-merge/lamp pipeline for TES detection algorithms
- minimap2 for splice-aware alignment
- scikit-learn for Gaussian Mixture Models

## Version History

- **v1.0.0** (2026-02-05): Initial release
  - GTF processing and transcript structure extraction
  - Minimap2 integration with junction guidance
  - GMM-based TES detection with consensus calling
  - Parallel processing support
