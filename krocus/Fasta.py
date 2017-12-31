'''Read in a FASTA file and extract all the k-mers'''
import operator
from Bio import SeqIO
from krocus.Kmers import Kmers

class Fasta:
	def __init__(self,logger, filename, k, num_top_kmers):
		self.logger = logger
		self.filename = filename
		self.k = k
		self.num_top_kmers = num_top_kmers
	
		self.sequences_to_kmers = self.sequence_kmers()
		self.all_kmers = self.all_kmers_in_file()
		self.top_kmers = self.top_kmers_in_file()
		self.avg_gene_length = 0

	def sequence_kmers(self):
		total_sequence_length = 0 
		seq_counter = 0
		
		kmer_to_sequences = {}
		for record in SeqIO.parse(self.filename, "fasta"):
			kmers = Kmers(str(record.seq), self.k)
			# We assume here that the sequence name is unique in the FASTA file
			kmer_to_sequences[record.id] = kmers.get_all_kmers()
			
			total_sequence_length  += len(record.seq)
			seq_counter += 1
			
		if seq_counter > 0:
			self.avg_gene_length  = total_sequence_length/seq_counter
			
		return kmer_to_sequences
		
	def all_kmers_in_file(self):
		all_kmers = {}
		for seq_name, kmer_counts in self.sequences_to_kmers.items():
			for kmer, count in kmer_counts.items():
				if kmer in all_kmers:
					all_kmers[kmer] += 1
				else:
					all_kmers[kmer] = 1
		return all_kmers
  	
	# With MLST most of the sequences will be near identical with just 1 SNP difference so keep a few top ones
	# as a quick filter
	def top_kmers_in_file(self):
		top = []
		
		loop_count = 0 
		for (k_key, k_count) in sorted(self.all_kmers.items(), key=operator.itemgetter(1), reverse=True):
			if loop_count >= self.num_top_kmers:
				break
			top.append(k_key)
			loop_count += 1 
		return top
		
		
		
		