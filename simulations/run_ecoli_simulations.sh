#!/bin/bash
# Iterate over all genomes and run simulations

# Install software
conda install -y -c bioconda nanosim-h pyfastaq ncbi-genome-download krocus mlst

# Download databases
ncbi-genome-download --format fasta --assembly-level complete --genus "Escherichia coli" bacteria
krocus_database_downloader -s "Escherichia coli#1" -o ecoli

declare -a errormodels=("ecoli_R9_1D" "ecoli_R9_2D")
READS=20000

cd refseq/bacteria
for i in $( ls -d * ); do 
	cd $i
	echo `pwd`
	
	# The MLST based on the assembled reference genome
	mlst --quiet *.fna.gz > ref.mlst
	cat ref.mlst
	TARGETST="$( awk '{print $3}' ref.mlst )"
	
	# simulate each error model
	for e in "${errormodels[@]}"; do
		# generate simulated reads in FASTA format
		gunzip -c  *.fna.gz | nanosim-h -n $READS -p ${e} - 
		
		# convert the FASTA to a FASTQ file
		fastaq to_fake_qual simulated_reads.fasta simulated_reads.fasta.qual
		fastaq fasta_to_fastq simulated_reads.fasta simulated_reads.fasta.qual simulated_reads.fastq
		
		# Calculate MLST from the FASTQ file
		{ time krocus --target_st ${TARGETST} -o ${e}.krocus ../../../ecoli simulated_reads.fastq ; } 2>${e}.time
		tail -n 1 ${e}.krocus
		
		# Cleanup
		rm simulated_reads.fasta simulated_reads.fastq simulated.log simulated_error_profile simulated_reads.fasta.qual

	done
	cd ..
done
