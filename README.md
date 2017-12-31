# Krocus
Given some raw uncorrected long reads, such as those from PacBio or Oxford Nanopore, calculate the multi-locus sequence type (MLST).

[![Build Status](https://travis-ci.org/andrewjpage/krocus.svg?branch=master)](https://travis-ci.org/andrewjpage/krocus)

# Installation
The only dependancy is Python3. Assuming you have python 3.3+ and pip installed, just run:
```
pip3 install git+git://github.com/andrewjpage/krocus.git
```

## Debian/Ubuntu (Trusty/Xenial)
To install Python3 on Ubuntu, as root run:
```
apt-get update -qq
apt-get install -y git python3 python3-setuptools python3-biopython python3-pip
pip3 install git+git://github.com/andrewjpage/krocus.git
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
                        Filename to save matching reads to
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        Output file [STDOUT]
  --min_fasta_hits MIN_FASTA_HITS, -m MIN_FASTA_HITS
                        Minimum No. of kmers matching a read
  --print_interval PRINT_INTERVAL, -p PRINT_INTERVAL
                        Print ST every this number of reads
  --kmer KMER, -k KMER  Kmer size [11]
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


