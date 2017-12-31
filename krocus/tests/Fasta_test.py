import unittest
import os
import logging
from krocus.Fasta import Fasta

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fasta')

class TestFasta(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		f = Fasta(logger, os.path.join(data_dir,'sample1.fa'),4, 1)
		self.assertEqual(f.sequence_kmers(),{'gene1': {'AAAA': 0}, 'gene2': {'TTTT': 0}, 'gene3': {'CCCC': 0}})
		self.assertEqual(f.all_kmers_in_file(), {'AAAA': 1, 'CCCC': 1, 'TTTT': 1})
		self.assertEqual(f.top_kmers, ['AAAA'])