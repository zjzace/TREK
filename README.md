# TREK - Transcription End site identifier

TREK identifies alternative polyadenylation (APA) sites and transcription end sites (TES) from long-read sequencing data.

## Installation

### Using Mamba/Conda

```bash
# Create environment
mamba env create -f environment.yml

# Activate environment
mamba activate TREK

# Install TREK
pip install -e .
```

### Using pip only

```bash
pip install -e .
```

## Requirements

- Python >= 3.9
- minimap2 (for alignment)
- Perl (for GTF to BED conversion)

## Usage

```bash
trek -g reference.gtf \
     -f genome.fa \
     -q reads.fastq.gz \
     -o results \
     -t 8
```

### Required Arguments

- `-g, --gtf`: Reference GTF annotation file
- `-f, --fasta`: Genome FASTA file  
- `-q, --fastq`: Input FASTQ files (can specify multiple)

### Optional Arguments

- `-o, --output`: Output directory (default: results)
- `-p, --prefix`: Output file prefix (default: polyA)
- `-t, --threads`: Number of threads (default: 8)
- `-j, --jobs`: Parallel jobs for TES detection (default: -1, all CPUs)

### TES Detection Parameters

- `--min-mapq`: Minimum mapping quality (default: 1)
- `--min-overlap`: Minimum overlap ratio (default: 0.2)
- `--min-reads`: Minimum reads per transcript (default: 10)
- `--min-cluster-size`: Minimum cluster size (default: 10)
- `--max-clusters`: Maximum clusters (default: 5)
- `--min-distance`: Minimum distance between sites in bp (default: 50)
- `--min-dominance`: Minimum relative dominance (default: 0.1)
- `--min-sharpness`: Minimum peak sharpness (default: 0.5)

### Internal Priming Filter Parameters

TREK automatically filters out APA sites that may result from internal priming artifacts. For transcripts with multiple APA sites, the filter examines the genomic sequence immediately upstream of each **non-dominant** site and removes those with high A-content. The dominant site (highest read count) is always kept.

- `--no-filter-priming`: Disable internal priming filter (enabled by default)
- `--priming-window`: Upstream window size (bp) for A-content check (default: 10)
- `--priming-a-threshold`: Maximum allowed A-content fraction (default: 0.5)

**How it works:**
1. For transcripts with multiple APA sites, the dominant site (highest abundance) is always kept
2. For non-dominant sites only, extract the upstream 10bp genomic sequence before each site
3. Calculate A proportion in the upstream window:
   - For + strand: counts A content in upstream sequence
   - For - strand: counts T content in downstream sequence (complement of A in transcript)
4. Remove non-dominant sites where A proportion > threshold (default 50%)
5. Recalculate abundances based on remaining sites

## Output Files

- `{prefix}.apa_sites.txt`: Main results with identified APA sites
- `{prefix}.summary.txt`: Summary statistics
- `{prefix}.read_assignments.pkl`: Read assignments (pickle format)
- `{prefix}.junctions.bed`: Junction guide for alignment

## Citation

If you use TREK in your research, please cite:

```
[Citation to be added]
```
