import unittest
import os
import logging
from krocus.Fastas import Fastas

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fastas')

class TestFastas(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		f = Fastas(logger, data_dir,4)
		self.assertEqual(list(f.fastas_to_kmers.values()), [{'AAAA': 1, 'TTTT': 1, 'CCCC': 1}, {'GGGG': 1, 'ATAT': 1, 'TATA': 1, 'CCCC': 1}])
