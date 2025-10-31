import unittest
import os
import logging
from krocus.Fastas import Fastas

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'fastas')

class TestFastas(unittest.TestCase):

    def test_allele_filenames(self):
        logger = logging.getLogger(__name__)
        f = Fastas(logger, data_dir, 4, False)
        filenames = f.filenames
        # Check that we found some .tfa files
        self.assertIsInstance(filenames, list)
        self.assertTrue(len(filenames) > 0)
        # All should end with .tfa
        for fn in filenames:
            self.assertTrue(fn.endswith('.tfa'))

    def test_fastas_to_kmers(self):
        logger = logging.getLogger(__name__)
        f = Fastas(logger, data_dir, 4, False)
        # Should have created a dictionary
        self.assertIsInstance(f.fastas_to_kmers, dict)
        # Each value should be a dictionary of kmers
        for fasta_obj, kmers in f.fastas_to_kmers.items():
            self.assertIsInstance(kmers, dict)

    def test_multiple_fasta_files(self):
        logger = logging.getLogger(__name__)
        f = Fastas(logger, data_dir, 4, False)
        # Should have processed multiple files
        self.assertGreater(len(f.fastas_to_kmers), 0)

    def test_divisible_by_3(self):
        logger = logging.getLogger(__name__)
        f = Fastas(logger, data_dir, 4, True)
        # Should still create the object even with filter enabled
        self.assertIsInstance(f.fastas_to_kmers, dict)

    def test_different_kmer_sizes(self):
        logger = logging.getLogger(__name__)
        for kmer_size in [5, 7, 9]:
            f = Fastas(logger, data_dir, kmer_size, False)
            self.assertIsInstance(f.fastas_to_kmers, dict)

if __name__ == '__main__':
    unittest.main()
