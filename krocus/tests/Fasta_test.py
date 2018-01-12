import unittest
import os
import logging
from krocus.Fasta import Fasta

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fasta')

class TestFasta(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		f = Fasta(logger, os.path.join(data_dir,'sample1.fa'),4, False)
		self.assertEqual(f.sequence_kmers(),{'gene1': {'GCAA': 0, 'CAAC': 0, 'AACG': 0, 'ACGC': 0, 'CGCT': 0, 'GCTG': 0, 'CTGA': 0, 'TGAC': 0, 'GACG': 0, 'ACGG': 0, 'CGGA': 0, 'GGAA': 0, 'GAAA': 0, 'AAAA': 0, 'AAAC': 0, 'ACGA': 0, 'CGAT': 0, 'GATC': 0, 'ATCT': 0, 'TCTG': 0, 'CTGG': 0, 'TGGT': 0, 'GGTT': 0, 'GTTT': 0, 'TTTT': 0, 'TTTG': 0, 'TTGC': 0, 'TGCC': 0, 'GCCC': 0, 'CCCT': 0, 'CCTT': 0, 'CTTT': 0, 'TTTC': 0, 'TTCA': 0, 'TCAC': 0, 'CACA': 0, 'ACAG': 0}, 'gene2': {'GGCC': 0, 'GCCG': 0, 'CCGC': 0, 'CGCA': 0, 'GCAC': 0, 'CACC': 0, 'ACCA': 0, 'CCAC': 0, 'CACG': 0, 'ACGG': 0, 'CGGC': 0, 'GGCG': 0, 'GCGC': 0, 'CGCT': 0, 'GCTC': 0, 'CTCG': 0, 'TCGC': 0, 'CGCC': 0, 'GCCC': 0, 'CCCT': 0, 'CCTT': 0, 'CTTC': 0}, 'gene3': {'CGCG': 0, 'GCGC': 0, 'CGCT': 0, 'GCTG': 0, 'CTGA': 0, 'TGAT': 0, 'GATT': 0, 'ATTT': 0, 'TTTT': 0, 'TTTG': 0, 'TTGC': 0, 'TGCG': 0, 'GCGT': 0, 'CGTG': 0, 'GTGG': 0, 'TGGC': 0, 'GGCA': 0, 'GCAA': 0, 'CAAT': 0, 'AATG': 0, 'ATGG': 0, 'GGCG': 0, 'GCGG': 0}})
		self.assertEqual(f.all_kmers_in_file(), {'GCAA': 2, 'CAAC': 1, 'AACG': 1, 'ACGC': 1, 'CGCT': 3, 'GCTG': 2, 'CTGA': 2, 'TGAC': 1, 'GACG': 1, 'ACGG': 2, 'CGGA': 1, 'GGAA': 1, 'GAAA': 1, 'AAAA': 1, 'AAAC': 1, 'ACGA': 1, 'CGAT': 1, 'GATC': 1, 'ATCT': 1, 'TCTG': 1, 'CTGG': 1, 'TGGT': 1, 'GGTT': 1, 'GTTT': 1, 'TTTT': 2, 'TTTG': 2, 'TTGC': 2, 'TGCC': 1, 'GCCC': 2, 'CCCT': 2, 'CCTT': 2, 'CTTT': 1, 'TTTC': 1, 'TTCA': 1, 'TCAC': 1, 'CACA': 1, 'ACAG': 1, 'GGCC': 1, 'GCCG': 1, 'CCGC': 1, 'CGCA': 1, 'GCAC': 1, 'CACC': 1, 'ACCA': 1, 'CCAC': 1, 'CACG': 1, 'CGGC': 1, 'GGCG': 2, 'GCGC': 2, 'GCTC': 1, 'CTCG': 1, 'TCGC': 1, 'CGCC': 1, 'CTTC': 1, 'CGCG': 1, 'TGAT': 1, 'GATT': 1, 'ATTT': 1, 'TGCG': 1, 'GCGT': 1, 'CGTG': 1, 'GTGG': 1, 'TGGC': 1, 'GGCA': 1, 'CAAT': 1, 'AATG': 1, 'ATGG': 1, 'GCGG': 1})
		
		
		
		
		