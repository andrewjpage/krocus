import unittest
import os
import logging
from krocus.Fastas import Fastas

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fastas')

class TestFastas(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		f = Fastas(logger, data_dir,4, False)
		self.assertEqual(list(f.fastas_to_kmers.values()), [{'CTGA': 1, 'TGAC': 1, 'GACG': 1, 'ACGC': 2, 'CGCC': 1, 'GCCA': 1, 'CCAA': 1, 'CAAG': 1, 'AAGG': 1, 'AGGT': 1, 'GGTG': 1, 'GTGG': 1, 'TGGC': 1, 'GGCG': 1, 'GCGG': 2, 'CGGA': 1, 'GGAG': 1, 'GAGG': 1, 'AGGC': 1, 'TTCG': 1, 'TCGT': 1, 'CGTC': 2, 'GTCG': 1, 'TCGC': 1, 'CGCG': 1, 'CGGC': 1, 'GGCT': 1, 'GCAA': 1, 'CAAC': 1, 'AACG': 1, 'CGCT': 1, 'GCTA': 1, 'CTAA': 1, 'TAAC': 1, 'ACGT': 1, 'GTCA': 1}, {'CTGA': 2, 'TGAC': 2, 'GACG': 1, 'ACGC': 1, 'CGCC': 1, 'GCCG': 1, 'CCGG': 1, 'CGGA': 1, 'GGAT': 1, 'GATG': 1, 'ATGC': 1, 'TGCT': 1, 'GCTG': 1, 'GACA': 1, 'ACAA': 1, 'CAAC': 1, 'AACT': 1, 'ACTG': 1, 'CTGG': 1, 'TGGC': 1, 'TTCG': 1, 'TCGT': 1, 'CGTC': 1, 'GTCG': 1, 'TCGC': 1, 'CGCT': 1, 'GCTT': 1, 'CTTC': 1, 'TTCT': 1, 'TCTG': 1, 'GACC': 1, 'ACCT': 1, 'CCTC': 1, 'CTCC': 1, 'TCCC': 1, 'CCCA': 1, 'CCAG': 1, 'CAGG': 1, 'AGGT': 1, 'GGTG': 1, 'GTGA': 1, 'TGAT': 1, 'GATC': 1, 'ATCC': 1, 'TCCT': 1, 'GCAA': 1, 'CAAT': 1, 'AATT': 1, 'ATTG': 1, 'TTGC': 1, 'TGCC': 1, 'GCCC': 1, 'CCCT': 1, 'CCTT': 1, 'CTTT': 1, 'TTTC': 1, 'TTCA': 1}])
