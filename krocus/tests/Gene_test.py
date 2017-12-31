import unittest
import os
from krocus.Gene import Gene

class TestGene(unittest.TestCase):

	def test_full_coverage(self):
		g = Gene('ABC_123', 10, 0)
		self.assertEqual(g.allele_number(), 123)
		self.assertEqual(g.allele_name(), 'ABC')
		self.assertTrue(g.is_full_coverage())
		self.assertEqual(str(g), 'ABC(123)')

	def test_no_coverage(self):
		g = Gene('ABC_123', 0, 10)
		self.assertFalse(g.is_full_coverage())
		self.assertEqual(str(g), 'ABC(123)*')
		
	def test_medium_coverage(self):
		g = Gene('ABC_123', 5, 5)
		self.assertFalse(g.is_full_coverage())
		self.assertEqual(str(g), 'ABC(123)*')
	
	def test_delimiters(self):
		for d in ['_','-','.']:
			g = Gene('ABC'+d+'123', 5 ,5)
			self.assertEqual(g.allele_number(), 123)
			self.assertEqual(g.allele_name(), 'ABC')
			self.assertEqual(str(g), 'ABC(123)*')
			
	def test_complex_gene_names(self):
		g = Gene('Csp_aspA.1', 1, 1)
		self.assertEqual(str(g), 'Csp_aspA(1)*')
		
		g = Gene('lipL41_3.2', 1, 1)
		self.assertEqual(str(g), 'lipL41_3(2)*')
		
		g = Gene('SH_1200.1', 1, 1)
		self.assertEqual(str(g), 'SH_1200(1)*')
