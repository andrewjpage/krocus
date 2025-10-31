import unittest
import os
import logging
import tempfile
from krocus.Fastq import Fastq
from krocus.Fastas import Fastas
from krocus.MlstProfile import MlstProfile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'fastq')

class TestFastq(unittest.TestCase):

    def test_fastq_initialization(self):
        """Test that Fastq object can be initialized with required parameters"""
        logger = logging.getLogger(__name__)
        profile_file = os.path.join(data_dir, 'profile.txt')
        
        if not os.path.exists(profile_file):
            self.skipTest("Profile file not found")
        
        mlst_profile = MlstProfile(profile_file, duplicate_warnings=False)
        fastas = Fastas(logger, data_dir, 11, False)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write('@read1\n')
            f.write('ACGTACGTACGTACGTACGTACGTACGT\n')
            f.write('+\n')
            f.write('IIIIIIIIIIIIIIIIIIIIIIIIIIII\n')
            fastq_file = f.name
        
        try:
            fastq = Fastq(logger, fastq_file, 11, fastas.get_fastas_to_kmers(), 
                         10, mlst_profile, 200, None, None)
            self.assertIsNotNone(fastq)
            self.assertEqual(fastq.k, 11)
            self.assertEqual(fastq.min_fasta_hits, 10)
        finally:
            if os.path.exists(fastq_file):
                os.unlink(fastq_file)

    def test_open_file_read_regular_file(self):
        """Test opening a regular FASTQ file"""
        logger = logging.getLogger(__name__)
        profile_file = os.path.join(data_dir, 'profile.txt')
        
        if not os.path.exists(profile_file):
            self.skipTest("Profile file not found")
            
        mlst_profile = MlstProfile(profile_file, duplicate_warnings=False)
        fastas = Fastas(logger, data_dir, 11, False)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write('@read1\n')
            f.write('ACGTACGTACGTACGTACGTACGTACGT\n')
            f.write('+\n')
            f.write('IIIIIIIIIIIIIIIIIIIIIIIIIIII\n')
            fastq_file = f.name
        
        try:
            fastq = Fastq(logger, fastq_file, 11, fastas.get_fastas_to_kmers(),
                         10, mlst_profile, 200, None, None)
            fh = fastq.open_file_read()
            self.assertIsNotNone(fh)
            fh.close()
        finally:
            if os.path.exists(fastq_file):
                os.unlink(fastq_file)

if __name__ == '__main__':
    unittest.main()
