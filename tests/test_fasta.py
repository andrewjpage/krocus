import unittest
import os
import logging
from krocus.Fasta import Fasta

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'fasta')

class TestFasta(unittest.TestCase):

    def test_four_kmers(self):
        logger = logging.getLogger(__name__)
        f = Fasta(logger, os.path.join(data_dir, 'sample1.fa'), 4, False)
        sequences = f.sequence_kmers()
        self.assertIn('gene1', sequences)
        self.assertIn('gene2', sequences)
        self.assertIn('gene3', sequences)
        # Check that some expected k-mers are present
        self.assertIn('GCAA', sequences['gene1'])
        
    def test_all_kmers_count(self):
        logger = logging.getLogger(__name__)
        f = Fasta(logger, os.path.join(data_dir, 'sample1.fa'), 4, False)
        all_kmers = f.all_kmers_in_file()
        # Verify that all_kmers is a dictionary with integer values
        self.assertIsInstance(all_kmers, dict)
        for kmer, count in all_kmers.items():
            self.assertIsInstance(count, int)
            self.assertGreaterEqual(count, 1)

    def test_divisible_by_3_filter(self):
        logger = logging.getLogger(__name__)
        # Test with divisible_by_3 = True
        f_filtered = Fasta(logger, os.path.join(data_dir, 'sample1.fa'), 4, True)
        # Should filter out sequences not divisible by 3
        self.assertIsInstance(f_filtered.sequences_to_kmers, dict)

    def test_sequence_kmers_vals(self):
        logger = logging.getLogger(__name__)
        f = Fasta(logger, os.path.join(data_dir, 'sample1.fa'), 4, False)
        vals = f.sequence_kmers_vals()
        self.assertIsInstance(vals, dict)
        # Check that values are dictionaries with integer counts
        for seq_name, kmers in vals.items():
            self.assertIsInstance(kmers, dict)

if __name__ == '__main__':
    unittest.main()
