# Krocus User Guide

## Overview

Krocus is a tool for rapid multi-locus sequence typing (MLST) directly from uncorrected long reads. It was designed to consume read data as it is produced, providing results in minutes. Krocus is particularly suited for long read sequencing technologies such as PacBio or Oxford Nanopore.

## Key Features

- **Fast**: Provides sequence types on average within 90 seconds
- **Direct from raw reads**: Works with uncorrected long reads (no assembly required)
- **Real-time analysis**: Can analyze data as it is being sequenced
- **High accuracy**: 94% sensitivity and 97% specificity
- **Multiple platforms**: Supports PacBio and Oxford Nanopore reads

## Installation

### Prerequisites

- Python 3.3 or higher
- pip (Python package installer)

### Install from PyPI

```bash
pip3 install krocus
```

### Install from source

```bash
pip3 install git+git://github.com/andrewjpage/krocus.git
```

### Platform-specific instructions

#### Debian/Ubuntu

```bash
apt-get update -qq
apt-get install -y git python3 python3-setuptools python3-biopython python3-pip
pip3 install krocus
```

#### Conda

```bash
conda install krocus
```

## Quick Start

### 1. Download MLST Database

First, download the MLST database for your species of interest:

```bash
# List available species
krocus_database_downloader -l

# Download database for a specific species
krocus_database_downloader --species "Salmonella enterica" --output_directory Salmonella_db
```

### 2. Run Krocus

Analyze your FASTQ file:

```bash
krocus Salmonella_db reads.fastq
```

## Detailed Usage

### krocus_database_downloader

This script downloads MLST databases from PubMLST.

```bash
krocus_database_downloader [options]
```

**Options:**

- `-l, --list_species`: List all available species
- `-s SPECIES, --species SPECIES`: Species to download
- `-o OUTPUT_DIR, --output_directory OUTPUT_DIR`: Output directory (default: mlst_files)
- `-v, --verbose`: Turn on debugging
- `--version`: Show version number

**Example:**

```bash
# List all available species
krocus_database_downloader --list_species

# Download Escherichia coli database
krocus_database_downloader --species "Escherichia coli" --output_directory ecoli_db
```

### krocus

The main analysis tool for MLST typing from reads.

```bash
krocus [options] allele_directory input.fastq
```

**Required Arguments:**

- `allele_directory`: Directory containing MLST database (from krocus_database_downloader)
- `input_fastq`: Input FASTQ file (can be gzipped, or use `-` for stdin)

**Key Options:**

- `-k KMER, --kmer KMER`: K-mer size (default: 11)
  - Critical parameter - should be set based on your read error rate
  - Rule of thumb: use the longest stretch of bases you can get without errors
  - For ~10% error rate, use k=9-11
  
- `-m MIN_HITS, --min_fasta_hits MIN_HITS`: Minimum k-mer matches (default: 10)
  - Lower values = more sensitive but slower
  - Higher values = faster but may miss some genes

- `-p INTERVAL, --print_interval INTERVAL`: Print results every N reads (default: 500)
  - Useful for real-time analysis
  - Lower values = more frequent updates

- `-f FILE, --filtered_reads_file FILE`: Save matching reads to file
  - Saves only regions predicted to contain MLST genes
  - Useful for downstream analysis

- `-o FILE, --output_file FILE`: Write results to file instead of stdout

- `-v, --verbose`: Enable debug output

**Advanced Options:**

- `--max_gap MAX_GAP`: Maximum gap for blocks (default: 4)
- `--margin MARGIN`: Flanking region around blocks (default: 50)
- `--min_block_size SIZE`: Minimum block size in bases (default: 150)
- `-d, --divisible_by_3`: Exclude genes not divisible by 3
- `--target_st ST`: For performance testing (print time to find given ST)

## Understanding K-mer Size

The k-mer size is the most important parameter. It directly relates to the error rate of your reads:

- **Higher k-mer size**: More specific, but requires lower error rates
- **Lower k-mer size**: More tolerant of errors, but less specific

**Guidelines:**
- Error rate ~5%: k = 15-19
- Error rate ~10%: k = 9-11
- Error rate ~15%: k = 7-9
- Error rate >20%: k = 5-7

## Real-time Analysis

Krocus can analyze data as it's being produced by the sequencer:

```bash
# Pipe data from a running sequencer
tail -f reads.fastq | krocus database_dir -

# Print results every 10 reads for faster updates
tail -f reads.fastq | krocus database_dir - -p 10
```

## Output Format

Krocus outputs tab-separated values:

```
ST    Coverage    Gene_Alleles
323   97.23       infB(1)  pgi(1)  phoE(9)*  tonB(93)  rpoB(1)*  gapA(2)  mdh(1)
```

**Fields:**

1. **ST**: Sequence type number
2. **Coverage**: Percentage of k-mers covered (0-100)
3. **Gene_Alleles**: Gene names with allele numbers in parentheses
   - `*` indicates partial match (some k-mers missing)

**Interpretation:**

- ST 323 identified
- 97.23% of expected k-mers found
- 2 genes have incomplete coverage (phoE and rpoB)
- Partial matches may indicate:
  - Read errors in that region
  - Insufficient read coverage
  - Novel allele

## Saved Reads

When using `--filtered_reads_file`, Krocus saves reads that match MLST genes:

```bash
krocus database_dir reads.fastq -f matching_reads.fastq
```

These filtered reads can be used for:
- De novo assembly
- Detailed analysis of specific genes
- Quality assessment
- Further bioinformatics analysis

## Performance Tips

1. **Adjust k-mer size** based on your read quality
2. **Use appropriate print interval** for your needs:
   - Real-time monitoring: `-p 10` or `-p 50`
   - Final result only: `-p 10000` (large number)
3. **Stream large files** instead of loading entirely:
   ```bash
   gunzip -c large_file.fastq.gz | krocus db_dir -
   ```

## Troubleshooting

### No results or low coverage

- **Check k-mer size**: May be too large for your error rate
- **Check database**: Ensure correct species database
- **Check read quality**: Very poor quality reads may not work
- **Try lower min_fasta_hits**: Default may be too stringent

### Results take too long

- **Increase k-mer size**: If reads are good quality
- **Increase min_fasta_hits**: More stringent filtering
- **Increase print_interval**: Reduce output overhead

### Partial matches (asterisks)

This is normal and expected, especially with:
- High error rate reads
- Low coverage regions
- Novel alleles

If >50% of genes are partial matches:
- Check k-mer size (may be too large)
- Check read coverage
- Verify correct database

## Resource Usage

- **Memory**: Approximately 1MB per 1MB of input FASTQ
  - Example: 550MB FASTQ → 550MB RAM
- **CPU**: Single-threaded processing
- **Speed**: ~90 seconds average for typical bacterial genome

## Citation

If you use Krocus in your research, please cite:

Andrew J. Page, Jacqueline A. Keane. (2018) Rapid multi-locus sequence typing direct from uncorrected long reads using Krocus. PeerJ 6:e5233 https://doi.org/10.7717/peerj.5233

Also consider citing PubMLST:

Keith A. Jolley, James E. Bray, Martin C. J. Maiden. (2018) Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Research. 3:124 https://doi.org/10.12688/wellcomeopenres.14826.1

## Support

- **Issues**: https://github.com/andrewjpage/krocus/issues
- **Email**: andrewjpage+krocus@gmail.com

## License

GNU GPL version 3
