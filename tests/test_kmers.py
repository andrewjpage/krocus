import unittest
import os
import logging
from krocus.Kmers import Kmers

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'kmers')

class TestKmers(unittest.TestCase):

    def test_four_kmers(self):
        k = Kmers('AAAAATTTTT', 4)
        self.assertEqual(k.get_all_kmers_counter(max_kmer_count=5), {'AAAA': 0, 'AAAT': 0, 'AATT': 0, 'ATTT': 0, 'TTTT': 0})
        
    def test_four_kmers_all(self):
        k = Kmers('AAAAATTTTTTTT', 4)
        self.assertEqual(k.get_all_kmers_freq(max_kmer_count=5), {'AAAA': 2, 'AAAT': 1, 'AATT': 1, 'ATTT': 1, 'TTTT': 4})
        
    def test_short_sequence(self):
        k = Kmers('A', 10)
        self.assertEqual(k.get_all_kmers_counter(), {})

    def test_empty_sequence(self):
        k = Kmers('', 4)
        self.assertEqual(k.get_all_kmers_counter(), {})

    def test_kmer_size_equal_to_sequence(self):
        k = Kmers('ATCG', 4)
        # When kmer size equals sequence length, no kmers are produced
        # because range(0, 0) is empty
        self.assertEqual(k.get_all_kmers_counter(max_kmer_count=5), {})

    def test_varying_kmer_sizes(self):
        for kmer_size in [5, 7, 9, 11]:
            k = Kmers('ATCGATCGATCGATCG', kmer_size)
            result = k.get_all_kmers_counter(max_kmer_count=10)
            self.assertIsInstance(result, dict)
            self.assertTrue(len(result) > 0)

    def test_get_all_kmers(self):
        k = Kmers('AAAAATTTTT', 4)
        result = k.get_all_kmers(max_kmer_count=5)
        self.assertIsInstance(result, dict)
        # Should contain KmerHit objects
        for kmer, hit in result.items():
            self.assertIsInstance(hit.count, int)
            self.assertIsInstance(hit.coordinates, list)

    def test_get_one_x_coverage_of_kmers(self):
        k = Kmers('ATCGATCGATCGATCG', 4)
        kmers = k.get_one_x_coverage_of_kmers()
        self.assertIsInstance(kmers, list)
        # Each kmer should be 4 bases long
        for kmer in kmers:
            self.assertEqual(len(kmer), 4)

if __name__ == '__main__':
    unittest.main()
