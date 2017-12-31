'''Read in a FASTQ file and identify matching alleles'''
from Bio import SeqIO
from krocus.Kmers import Kmers
from krocus.Read import Read
from krocus.Gene import Gene
from krocus.Blocks import Blocks
import subprocess
import numpy

class Error (Exception): pass

class Fastq:
	def __init__(self,logger, filename, k, fasta_kmers, min_fasta_hits, mlst_profile, print_interval, output_file, filtered_reads_file):
		self.logger = logger
		self.filename = filename
		self.k = k
		self.fasta_kmers = fasta_kmers
		self.min_fasta_hits = min_fasta_hits
		self.mlst_profile = mlst_profile
		self.print_interval = print_interval
		self.output_file = output_file
		self.filtered_reads_file = filtered_reads_file
		self.max_gap = 4 # multiples of the kmer
		self.min_block_size = 150
		self.margin = 100

	def initial_read_filter(self):
		counter = 0 
		match_counter =0
		
		fh = self.open_file_read()
		read = Read()
		
		while read.get_next_from_file(fh):
			counter += 1
			self.map_kmers_to_read(read.seq, read)
			if counter % self.print_interval == 0:
				self.full_gene_coverage()
				
		self.full_gene_coverage()
								
		return self
		
	def map_kmers_to_read(self, sequence,read):
		seq_length = len(sequence)
		end = seq_length - self.k
		
		kmers_obj = Kmers(sequence, self.k)
		read_kmers = kmers_obj.get_all_kmers_array()
		is_read_matching = False
		
		for (fasta_obj, fasta_kmers) in self.fasta_kmers.items():
			
			sequence_hits = numpy.zeros(int(seq_length/self.k)+1, dtype=int)
			hit_counter = 0
			for i in range(0, end):
				if read_kmers[i] in fasta_kmers:
					hit_counter += 1
					sequence_hits[int(i/self.k)] += 1
			
			if hit_counter < self.min_fasta_hits:
				continue
				
			blocks_obj = Blocks(self.k, self.min_block_size, self.max_gap, self.margin)
			block_start, block_end = blocks_obj.find_largest_block(sequence_hits)
			block_start = blocks_obj.adjust_block_start(block_start)

			if block_end == 0:
				continue
			block_end = blocks_obj.adjust_block_end(block_end,seq_length)

			block_kmers = self.create_kmers_for_block(block_start, block_end, sequence)
			self.apply_kmers_to_genes(fasta_obj,block_kmers)
			is_read_matching = True
			
			if self.filtered_reads_file:
				with open(self.filtered_reads_file, 'a+') as output_fh:
					output_fh.write(str(read.subsequence(block_start, block_end)))
			
			
	def create_kmers_for_block(self, block_start, block_end, sequence):
		if block_end ==  0:
			return {}
		
		block_seq = sequence[block_start:block_end]
			
		kmers_obj = Kmers(block_seq, self.k)
		return kmers_obj.get_all_kmers()
	
		# {'gene1': {'AAAA': 1}, 'gene2': {'TTTT': 1}, 'gene3': {'CCCC': 1}}
	def apply_kmers_to_genes(self,fasta_obj,hit_kmers):
		for (gene_name, kmers_dict) in fasta_obj.sequences_to_kmers.items():
			for kmer in kmers_dict.keys():
				if kmer in hit_kmers:
					fasta_obj.sequences_to_kmers[gene_name][kmer] += 1 
		
	def full_gene_coverage(self):
		alleles = []
		for fasta_obj in self.fasta_kmers.keys():
			largest_gene = 0
			largest_gene_name = ''
			largest_zero = 0
			for (gene_name, kmers_dict) in fasta_obj.sequences_to_kmers.items():
				kv = kmers_dict.values()
				
				kv_coverage = [x for x in kv if x >= 1]
				kv_zero = [x for x in kv if x == 0]
				kl  = len(kv_coverage)
				kz = len(kv_zero)
				if kl > largest_gene:
					largest_gene = kl
					largest_gene_name = gene_name
					largest_zero = kz
	
			if largest_gene_name != '':
				alleles.append(Gene(largest_gene_name, largest_gene, largest_zero))
			
		gene_to_allele_number = {}
		for a in alleles:
			gene_to_allele_number[a.allele_name()] = a.allele_number()
		st = self.mlst_profile.get_sequence_type(gene_to_allele_number)
		self.output_st_and_alleles(st, alleles)
		
	def output_st_and_alleles(self, st, alleles):
		allele_string = str(st)+"\t" 
		for a in alleles:
			allele_string += str(a)+"\t"		
		if self.output_file:
			with open(self.output_file, 'a+') as output_fh:
				output_fh.write(allele_string + "\n")			
		else:
			print(allele_string)	
		# minority variants
	
	# Derived from https://github.com/sanger-pathogens/Fastaq
	# Author: Martin Hunt	
	def open_file_read(self):
		if self.filename == '-':
			f = sys.stdin
		elif self.filename.endswith('.gz'):
			# first check that the file is OK according to gunzip
			retcode = subprocess.call('gunzip -t ' + self.filename, shell=True)
			
			# now open the file
			f = os.popen('gunzip -c ' + self.filename)
		else:
			try:
				f = open(self.filename)
			except:
				raise Error("Error opening for reading file '" + self.filename + "'")
		
		return f
	
