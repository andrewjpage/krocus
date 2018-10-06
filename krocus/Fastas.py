'''Read multiple FASTA files and extract all the k-mers'''
import os
from Bio import SeqIO
from krocus.Fasta import Fasta

class Fastas:
	def __init__(self,logger, allele_directory, k, divisible_by_3, max_kmers = 10):
		self.logger = logger
		self.allele_directory = allele_directory
		self.divisible_by_3 = divisible_by_3
		self.filenames = self.allele_filenames(allele_directory)
		self.k = k
		self.max_kmers = max_kmers
		self.fastas_to_kmers = self.get_fastas_to_kmers()

	def get_fastas_to_kmers(self):
		fastas_to_kmers = {}
		for f in self.filenames:
			fasta_obj = Fasta(self.logger, f, self.k,self.divisible_by_3, max_kmers = self.max_kmers)
			fastas_to_kmers[fasta_obj] = fasta_obj.all_kmers
		return fastas_to_kmers

	def allele_filenames(self,allele_directory):
		fasta_filenames = []
		for dirpath, dirnames, filenames in os.walk(allele_directory):
		    for filename in [f for f in filenames if f.endswith(".tfa")]:
		        fasta_filenames.append(os.path.join(dirpath, filename))
		return fasta_filenames
		