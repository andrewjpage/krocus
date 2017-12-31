import logging
import os
from krocus.Fastas import Fastas
from krocus.Fastq import Fastq
from krocus.MlstProfile import MlstProfile

class Krocus:
	def __init__(self,options):
		self.logger = logging.getLogger(__name__)
		self.allele_directory           = options.allele_directory 
		self.input_fastq                = options.input_fastq
		self.kmer                       = options.kmer
		self.verbose                    = options.verbose
		self.min_fasta_hits             = options.min_fasta_hits
		self.print_interval             = options.print_interval
		self.num_top_kmers              = 5 # may not be needed
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
	def run(self):
		mlst_profile = MlstProfile(self.mlst_profile_file())
		fastas = Fastas(self.logger, self.allele_directory, self.kmer, self.num_top_kmers)
		fastq = Fastq(self.logger, self.input_fastq, self.kmer, fastas.get_fastas_to_kmers(), self.min_fasta_hits , mlst_profile, self.print_interval)
		fastq.initial_read_filter()

	def mlst_profile_file(self):
		profile_txt = os.path.join(self.allele_directory, 'profile.txt')
		if not os.path.exists(profile_text):
			self.logger.error("The MLST profile file cannot be accessed: "+ profile_txt)
		return profile_txt