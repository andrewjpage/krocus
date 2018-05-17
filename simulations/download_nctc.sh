#!/bin/bash
# Raw data accession number such as ERR1466818
ERRACCESSION=$1

# ensure tools are installed
conda install -c bioconda pbh5tools krocus mlst

# Download the raw data
wget http://sra-download.ncbi.nlm.nih.gov/srapub_files/${ERRACCESSION}_${ERRACCESSION}_hdf5.tgz

# unzip
tar -xvfz ${ERRACCESSION}_${ERRACCESSION}_hdf5.tgz

# convert 
bash5tools.py --outFilePrefix ${ERRACCESSION} --outType fastq --readType Raw *.bas.h5 

