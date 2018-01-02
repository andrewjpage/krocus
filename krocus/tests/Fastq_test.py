import unittest
import os
import logging
import filecmp
from krocus.Fastq import Fastq
from krocus.Fastas import Fastas
from krocus.MlstProfile import MlstProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fastq')

import cProfile, pstats, io

class TestFastq(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		fastas = Fastas(logger, data_dir,4, False)
		mlst_profile = MlstProfile(data_dir+'/profile.txt')
		
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 4 , fastas.get_fastas_to_kmers(), 1,  mlst_profile, 100, None, None)
				
		self.assertTrue(fastq.initial_read_filter())
		
	def test_salmonella_welt_pacbio(self):
		logger = logging.getLogger(__name__)
		kmer = 11
		fastas = Fastas(logger, data_dir,kmer, False)
		mlst_profile = MlstProfile(data_dir+'/profile.txt')

		pr = cProfile.Profile()
		pr.enable()
		
		fastq = Fastq(logger, os.path.join(data_dir,'pacbio.fastq'), kmer , fastas.get_fastas_to_kmers(), 5,  mlst_profile, 100, None, None)
		fastq.initial_read_filter()
		
		pr.disable()
		s = io.StringIO()
		sortby = 'cumulative'
		ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
		ps.print_stats()
		print(s.getvalue())
		# ST	aroC	dnaN	hemD	hisD	purE	sucA	thrA
		# 365	130	    97	    25	    125   	84	    9	    101


		
		
		
		