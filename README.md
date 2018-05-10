# Krocus
Genome sequencing is rapidly being adopted in reference labs and hospitals for bacterial outbreak investigation and diagnostics where time is critical. Seven gene multi-locus sequence typing is a standard tool for broadly classifying samples into sequence types, allowing, in many cases, to rule a sample in or out of an outbreak, or allowing for general characteristics about a bacterial strain to be inferred. Long read sequencing technologies, such as from PacBio or Oxford Nanopore, can produce read data within minutes of an experiment starting, unlike short read sequencing technologies which require many hours/days. However, the error rates of raw uncorrected long read data are very high. We present Krocus which can predict a sequence type directly from uncorrected long reads, and which was designed to consume read data as it is produced, providing results in minutes. It is the only tool which can do this from uncorrected long reads. We tested Krocus on over 600 samples sequenced with using long read sequencing technologies from PacBio and Oxford Nanopore. It provides sequence types on average within 90 seconds, with a sensitivity of 94% and specificity of 97%, directly from uncorrected raw sequence reads. The software is written in Python and is available under the open source license GNU GPL version 3.

[![Build Status](https://travis-ci.org/andrewjpage/krocus.svg?branch=master)](https://travis-ci.org/andrewjpage/krocus)

# Paper
"Rapid multi-locus sequence typing direct from uncorrected long reads using Krocus", Andrew J. Page, Jacqueline A. Keane, bioRxiv 259150; (2018) doi: https://doi.org/10.1101/259150

# Installation
The only dependancy is Python3. Assuming you have python 3.3+ and pip installed, just run:
```
pip3 install krocus
```

or if you wish to install the latest development version:
```
pip3 install git+git://github.com/andrewjpage/krocus.git
```

On some systems pip3 may be just called pip.

## Debian/Ubuntu (Trusty/Xenial)
To install Python3 on Ubuntu, as root run:
```
apt-get update -qq
apt-get install -y git python3 python3-setuptools python3-biopython python3-pip
pip3 install krocus
```

## Conda
First install Conda (Python3), then run:
```
conda install  krocus
```

## Windows
Like virtually all Bioinformatics software, this software is unlikely to work on Windows. Try using a Linux virtual machine.

# Usage
## krocus_database_downloader script
First of all you need MLST databases. There is a snapshot bundled with this repository for your convenience, or alternatively you can use the downloader script to get the latest data. You will need internet access for this step.

```
usage: krocus_database_downloader [options]

Download

optional arguments:
  -h, --help            show this help message and exit
  --list_species, -l    List all available species (default: False)
  --species SPECIES, -s SPECIES
                        Species to download (default: None)
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Output directory (default: mlst_files)
  --verbose, -v         Turn on debugging (default: False)
  --version             show program's version number and exit

```
First of all you can get a list of available databases by running:
```
krocus_database_downloader -l
```

From this list choose one of the species and use it for the next step:
```
krocus_database_downloader  --species "Salmonella enterica" --output_directory Salmonella_enterica
```
You will now have a directory called __Salmonella_enterica___ which can be provided to the main script.

## krocus script
This is the main script of the application. The manditory inputs are a directory containing an MLST database (from the previous step), and a FASTQ file, which can be optionally gzipped.
```
usage: krocus [options] allele_directory input.fastq

multi-locus sequence typing (MLST) from uncorrected long reads

positional arguments:
  allele_directory      Allele directory
  input_fastq           Input FASTQ file (optionally gzipped)

optional arguments:
  -h, --help            show this help message and exit
  --filtered_reads_file FILTERED_READS_FILE, -f FILTERED_READS_FILE
                        Filename to save matching reads to (default: None)
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        Output file [STDOUT] (default: None)
  --max_gap MAX_GAP     Maximum gap for blocks to be contigous, measured in
                        multiples of the k-mer size (default: 4)
  --margin MARGIN       Flanking region around a block to use for mapping
                        (default: 100)
  --min_block_size MIN_BLOCK_SIZE
                        Minimum block size in bases (default: 150)
  --min_fasta_hits MIN_FASTA_HITS, -m MIN_FASTA_HITS
                        Minimum No. of kmers matching a read (default: 10)
  --print_interval PRINT_INTERVAL, -p PRINT_INTERVAL
                        Print ST every this number of reads (default: 200)
  --kmer KMER, -k KMER  k-mer size (default: 11)
  --target_st TARGET_ST
                        For performance testing print time to find given ST
                        (default: None)
  --verbose, -v         Turn on debugging [0]
  --version             show program's version number and exit
```

### Required
__allele_directory__: The directory containing the MLST database you wish to query against. This is generated by the krocus_database_downloader script and just contains copies of the allele sequences in FASTA format and the profile.txt file linking allele numbers to STs.

__input_fastq__: This is a single FASTQ file. It can be optionally gzipped. Alternatively input can be read from stdin by using the dash character (-) as the input file name.

### Options
__kmer__:  The most important parameter. Long reads have a high error rate, so if you set this too high, nothing will match (because it will contain errors). If you set it too low, everything will match, which isnt much use to you. Thinking about your data, on average how long of a stretch of bases can you get in your read without errors? This is what you should set your kmer to. For example, if you have an average of 1 error every 10 bases, then the ideal kmer would be 9.

__min_fasta_hits__: This is the minimum number of matching kmers in a read, for the read to be considered for analysis. It is a hard minimum threshold which is really there to speed things up. If you set this too high, then nothing will be returned.

__filtered_reads_file__: If you provide a filename for this option, all of the reads which are estimated to match one of the MLST genes are saved to a file. Only the region predicted to contain the MLST gene is saved. This can be used for downstream analysis, such as de novo assembly. This file should not already exist. 

__output_file__: By default the predicted sequence types are printed to screen (STDOUT). If a filename is provided, the predicted sequence types are instead printed to this file.  This file should not already exist. 

__print_interval__: Print out the predicted sequence type every X number of reads. This is where you are performing analysis in real time and want a quick result.

# Output
The output format is: the predicted sequence type (ST) number (column 1), the percentage k-mer coverage of the alleles (0-100) (column 2), the specific alleles identified. For each allele the name of the gene is noted and the allele number, which corresponds to a unique sequence, is given in brackets. If there is only a partial match a start (*) is appended.

```
323	97.23	infB(1)	pgi(1)	phoE(9)*	tonB(93)	rpoB(1)*	gapA(2)	mdh(1)
```
In the above example the sequence type is 323,and 97.23% of k-mers making up 323 are covered by reads. 2 of the alleles are have k-mers which were not identified in the reads, possibly due to errors in the reads encountered or no reads covering that region were found. 



# Resource usage
For an 550Mbyte FASTQ file (unzipped) of long reads from a Pacbio RSII containing Salmonella required 550MB of RAM.


