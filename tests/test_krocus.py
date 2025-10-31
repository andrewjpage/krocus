import unittest
import os
import tempfile
import logging
from unittest.mock import Mock, patch
from krocus.Krocus import Krocus

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data')

class TestKrocus(unittest.TestCase):

    def test_initialization(self):
        """Test Krocus object initialization"""
        options = Mock()
        options.allele_directory = data_dir
        options.input_fastq = '/tmp/test.fastq'
        options.kmer = 11
        options.verbose = False
        options.min_fasta_hits = 10
        options.print_interval = 200
        options.output_file = None
        options.filtered_reads_file = None
        options.target_st = None
        options.max_gap = 4
        options.min_block_size = 150
        options.margin = 50
        options.divisible_by_3 = False
        options.min_kmers_for_onex_pass = 10
        options.max_kmers = 10
        
        krocus = Krocus(options)
        self.assertEqual(krocus.kmer, 11)
        self.assertEqual(krocus.min_fasta_hits, 10)
        self.assertEqual(krocus.allele_directory, data_dir)

    def test_output_file_exists_error(self):
        """Test that error is raised if output file already exists"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            output_file = f.name
        
        try:
            options = Mock()
            options.allele_directory = data_dir
            options.input_fastq = '/tmp/test.fastq'
            options.kmer = 11
            options.verbose = False
            options.min_fasta_hits = 10
            options.print_interval = 200
            options.output_file = output_file
            options.filtered_reads_file = None
            options.target_st = None
            options.max_gap = 4
            options.min_block_size = 150
            options.margin = 50
            options.divisible_by_3 = False
            options.min_kmers_for_onex_pass = 10
            options.max_kmers = 10
            
            with self.assertRaises(SystemExit):
                krocus = Krocus(options)
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_filtered_reads_file_exists_error(self):
        """Test that error is raised if filtered reads file already exists"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            filtered_file = f.name
        
        try:
            options = Mock()
            options.allele_directory = data_dir
            options.input_fastq = '/tmp/test.fastq'
            options.kmer = 11
            options.verbose = False
            options.min_fasta_hits = 10
            options.print_interval = 200
            options.output_file = None
            options.filtered_reads_file = filtered_file
            options.target_st = None
            options.max_gap = 4
            options.min_block_size = 150
            options.margin = 50
            options.divisible_by_3 = False
            options.min_kmers_for_onex_pass = 10
            options.max_kmers = 10
            
            with self.assertRaises(SystemExit):
                krocus = Krocus(options)
        finally:
            if os.path.exists(filtered_file):
                os.unlink(filtered_file)

    def test_mlst_profile_file(self):
        """Test mlst_profile_file method"""
        # Create a temporary directory with profile.txt
        temp_dir = tempfile.mkdtemp()
        profile_path = os.path.join(temp_dir, 'profile.txt')
        
        try:
            with open(profile_path, 'w') as f:
                f.write('ST\tabc\tdef\n')
                f.write('1\t1\t1\n')
            
            options = Mock()
            options.allele_directory = temp_dir
            options.input_fastq = '/tmp/test.fastq'
            options.kmer = 11
            options.verbose = False
            options.min_fasta_hits = 10
            options.print_interval = 200
            options.output_file = None
            options.filtered_reads_file = None
            options.target_st = None
            options.max_gap = 4
            options.min_block_size = 150
            options.margin = 50
            options.divisible_by_3 = False
            options.min_kmers_for_onex_pass = 10
            options.max_kmers = 10
            
            krocus = Krocus(options)
            profile_file = krocus.mlst_profile_file()
            self.assertEqual(profile_file, profile_path)
        finally:
            if os.path.exists(profile_path):
                os.unlink(profile_path)
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_verbose_logging(self):
        """Test that verbose mode sets DEBUG logging level"""
        options = Mock()
        options.allele_directory = data_dir
        options.input_fastq = '/tmp/test.fastq'
        options.kmer = 11
        options.verbose = True
        options.min_fasta_hits = 10
        options.print_interval = 200
        options.output_file = None
        options.filtered_reads_file = None
        options.target_st = None
        options.max_gap = 4
        options.min_block_size = 150
        options.margin = 50
        options.divisible_by_3 = False
        options.min_kmers_for_onex_pass = 10
        options.max_kmers = 10
        
        krocus = Krocus(options)
        self.assertEqual(krocus.logger.level, logging.DEBUG)

if __name__ == '__main__':
    unittest.main()
