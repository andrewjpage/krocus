import unittest
import argparse
import os
import tempfile
from krocus.InputTypes import InputTypes

class TestInputTypes(unittest.TestCase):

    def test_is_fastq_file_valid_existing_file(self):
        """Test validation of existing FASTQ file"""
        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write('@read1\nACGT\n+\nIIII\n')
            temp_file = f.name
        
        try:
            result = InputTypes.is_fastq_file_valid(temp_file)
            self.assertEqual(result, temp_file)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)

    def test_is_fastq_file_valid_stdin(self):
        """Test validation of stdin input"""
        result = InputTypes.is_fastq_file_valid('-')
        self.assertEqual(result, '-')

    def test_is_fastq_file_valid_nonexistent(self):
        """Test that error is raised for non-existent file"""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_fastq_file_valid('/nonexistent/file.fastq')

    def test_is_allele_directory_valid_existing(self):
        """Test validation of existing directory"""
        # Create a temporary directory
        temp_dir = tempfile.mkdtemp()
        
        try:
            result = InputTypes.is_allele_directory_valid(temp_dir)
            self.assertEqual(result, temp_dir)
        finally:
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_is_allele_directory_valid_nonexistent(self):
        """Test that error is raised for non-existent directory"""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_allele_directory_valid('/nonexistent/directory')

    def test_is_kmer_valid_valid_values(self):
        """Test valid k-mer values"""
        for kmer in [5, 7, 9, 11, 15, 20, 25, 31]:
            result = InputTypes.is_kmer_valid(str(kmer))
            self.assertEqual(result, kmer)

    def test_is_kmer_valid_boundary_values(self):
        """Test boundary k-mer values"""
        # Minimum valid value
        self.assertEqual(InputTypes.is_kmer_valid('5'), 5)
        # Maximum valid value
        self.assertEqual(InputTypes.is_kmer_valid('31'), 31)

    def test_is_kmer_valid_invalid_too_small(self):
        """Test that error is raised for k-mer < 5"""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid('4')
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid('0')

    def test_is_kmer_valid_invalid_too_large(self):
        """Test that error is raised for k-mer > 31"""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid('32')
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid('100')

    def test_is_kmer_valid_invalid_non_numeric(self):
        """Test that error is raised for non-numeric values"""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid('abc')
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid('10.5')

if __name__ == '__main__':
    unittest.main()
