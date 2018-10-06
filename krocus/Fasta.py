'''Read in a FASTA file and extract all the k-mers'''
import operator
from Bio import SeqIO
from krocus.Kmers import Kmers

class Fasta:
	def __init__(self,logger, filename, k, divisible_by_3, max_kmers = 5):
		self.logger = logger
		self.filename = filename
		self.max_kmers = max_kmers
		self.k = k
		self.divisible_by_3 = divisible_by_3
	
		self.sequences_to_kmers = self.sequence_kmers()
		self.sequences_to_kmers_count = self.sequence_kmers_vals()
		self.all_kmers = self.all_kmers_in_file()
		
		
	def sequence_kmers_vals(self):
		seq_counter = 0
	
		kmer_to_sequences = {}
		for record in SeqIO.parse(self.filename, "fasta"):
			sequence_length  = len(record.seq)
			if self.divisible_by_3 and sequence_length % 3 != 0:
				self.logger.warning("Excluding gene as it is not divisible by 3:"+record.id)
				continue
		
			kmers = Kmers(str(record.seq), self.k)
			# We assume here that the sequence name is unique in the FASTA file
			kmer_to_sequences[record.id] = kmers.get_all_kmers_freq(max_kmer_count = self.max_kmers)
		
			seq_counter += 1
		
		return kmer_to_sequences

	def sequence_kmers(self):
		seq_counter = 0
		
		kmer_to_sequences = {}
		for record in SeqIO.parse(self.filename, "fasta"):
			sequence_length  = len(record.seq)
			if self.divisible_by_3 and sequence_length % 3 != 0:
				self.logger.warning("Excluding gene as it is not divisible by 3:"+record.id)
				continue
			
			kmers = Kmers(str(record.seq), self.k)
			# We assume here that the sequence name is unique in the FASTA file
			kmer_to_sequences[record.id] = kmers.get_all_kmers_counter(max_kmer_count = self.max_kmers)
			
			seq_counter += 1
			
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

		
		
		
		