import unittest
import os
import logging
import filecmp
from krocus.Fastq import Fastq
from krocus.Fastas import Fastas
from krocus.MlstProfile import MlstProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fastq')
databases =  os.path.join(test_modules_dir, '..','..','databases')

import cProfile, pstats, io

class TestFastq(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		fastas = Fastas(logger, data_dir,4, False)
		mlst_profile = MlstProfile(data_dir+'/profile.txt')
		
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 4 , fastas.get_fastas_to_kmers(), 1,  mlst_profile, 50, None, None)
				
		self.assertTrue(fastq.read_filter_and_map())
		
	def test_salmonella_welt_pacbio(self):
		logger = logging.getLogger(__name__)
		kmer = 11
		fastas = Fastas(logger, data_dir,kmer, False)
		mlst_profile = MlstProfile(data_dir+'/profile.txt')

		fastq = Fastq(logger, os.path.join(data_dir,'pacbio.fastq.gz'), kmer , fastas.get_fastas_to_kmers(), 10,  mlst_profile, 50, None, None)
		self.assertTrue(fastq.read_filter_and_map())

	def test_stap_aureus_pacbio(self):
		logger = logging.getLogger(__name__)
		kmer = 11
		
		pr = cProfile.Profile()
		pr.enable()      
		
		fastas = Fastas(logger, databases+'/Staphylococcus_aureus/',kmer, False)
		mlst_profile = MlstProfile(databases+'/Staphylococcus_aureus/profile.txt')
		
		fastq = Fastq(logger, os.path.join(data_dir,'NCTC11150.fastq.gz'), kmer , fastas.get_fastas_to_kmers(), 10,  mlst_profile, 50, None, None, max_kmers = 1)
		self.assertTrue(fastq.read_filter_and_map())

		pr.disable()
		s = io.StringIO()
		sortby = 'cumulative'
		ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
		ps.print_stats()
		print(s.getvalue())

	def test_kp_nanopore(self):
		logger = logging.getLogger(__name__)
		kmer = 11
		fastas = Fastas(logger, databases+'/Klebsiella_pneumoniae/',kmer, False)
		mlst_profile = MlstProfile(databases+'/Klebsiella_pneumoniae/profile.txt')
		
		fastq = Fastq(logger, os.path.join(data_dir,'nanopore_kp_12.fastq.gz'), kmer , fastas.get_fastas_to_kmers(), 10,  mlst_profile, 50, None, None)
		self.assertTrue(fastq.read_filter_and_map())

