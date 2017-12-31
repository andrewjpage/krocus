'''Read in a FASTQ file and identify matching alleles'''
from Bio import SeqIO
from krocus.Kmers import Kmers
from krocus.Read import Read
from krocus.Gene import Gene
import subprocess
import numpy

class Error (Exception): pass

class Fastq:
	def __init__(self,logger, filename, k, fasta_kmers, min_fasta_hits, mlst_profile, print_interval):
		self.logger = logger
		self.filename = filename
		self.k = k
		self.fasta_kmers = fasta_kmers
		self.min_fasta_hits = min_fasta_hits
		self.mlst_profile = mlst_profile
		self.print_interval = print_interval
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
			self.map_kmers_to_read(read.seq)
			if counter % self.print_interval == 0:
				self.full_gene_coverage()
								
		return self
		
		
	def map_kmers_to_read(self, sequence):
		seq_length = len(sequence)
		end = seq_length - self.k
		
		kmers_obj = Kmers(sequence, self.k)
		read_kmers = kmers_obj.get_all_kmers_array()
		
		for (fasta_obj, fasta_kmers) in self.fasta_kmers.items():
			
			sequence_hits = numpy.zeros(int(seq_length/self.k)+1, dtype=int)
			hit_counter = 0
			for i in range(0,end):
				if read_kmers[i] in fasta_kmers:
					hit_counter += 1
					sequence_hits[int(i/self.k)] += 1
			
			if hit_counter < self.min_fasta_hits:
				continue
				
			block_start, block_end = self.find_largest_block(sequence_hits)
			block_start *= self.k
			block_end *= self.k
			if block_start - self.margin < 0:
				block_start = 0
			else:
				block_start -= self.margin
			
			if block_end == 0:
				continue
			elif block_end + self.margin > seq_length:
				block_end = seq_length
			else:
				block_end +=  self.margin

			block_kmers = self.create_kmers_for_block(block_start, block_end, sequence)
			self.apply_kmers_to_genes(fasta_obj,block_kmers)
			
				
	def create_kmers_for_block(self, block_start, block_end, sequence):
		if block_end ==  0:
			return {}
		
		block_seq = sequence[block_start:block_end]
			
		kmers_obj = Kmers(block_seq, self.k)
		return kmers_obj.get_all_kmers()
				
		
	def merge_blocks(self, blocks):
		for i in range(0,len(blocks)-1 ):
			if blocks[i][1] + self.max_gap > blocks[i+1][0]:
				if blocks[i][1]  < blocks[i+1][1]:
					blocks[i][1] = 	blocks[i+1][1]
				blocks[i+1][0] = blocks[i][0]
				blocks[i+1][1] = blocks[i][1]
		return blocks
		
	
	def find_largest_block(self, sequence_hits):
		blocks = self.find_all_blocks(sequence_hits)
		merged_blocks = self.merge_blocks(blocks)
		
		largest_block = 0
		largest_block_index = 0
		
		for i, block in enumerate(merged_blocks):
			block_size = block[1] - block[0]
			if block_size > largest_block:
				largest_block_index = i
				largest_block = block_size
				
		if largest_block < (self.min_block_size/self.k):
			return 0,0
				
		return merged_blocks[largest_block_index][0], merged_blocks[largest_block_index][1]
		
	def find_all_blocks(self, sequence_hits):
		blocks = []
		in_block = False
		current_block_start = 0
		for i,val_count in enumerate(sequence_hits):
			
			if not in_block and val_count > 0:
				in_block = True
				current_block_start = i
			elif in_block and val_count == 0: 
				in_block = False
				blocks.append([current_block_start, i])
				
		if in_block:
			blocks.append([current_block_start, len(sequence_hits)])
				
		return blocks	
	
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
		print(str(st)+"\t", end='')
		for a in alleles:
			print(str(a)+"\t", end='')
		print("")
	
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
	
