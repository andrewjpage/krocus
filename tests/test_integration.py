"""
Integration tests for Krocus
These tests verify the interaction between different modules.
"""
import unittest
import os
import tempfile
import logging
from krocus.Fastas import Fastas
from krocus.MlstProfile import MlstProfile
from krocus.Fasta import Fasta
from krocus.Kmers import Kmers
from krocus.Gene import Gene

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data')

class TestIntegration(unittest.TestCase):

    def test_fastas_and_mlst_profile_integration(self):
        """Test that Fastas and MlstProfile work together"""
        logger = logging.getLogger(__name__)
        
        # Create a temporary profile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\n')
            f.write('1\t1\t1\n')
            f.write('2\t1\t2\n')
            profile_file = f.name
        
        # Create temporary fasta directory
        temp_dir = tempfile.mkdtemp()
        fasta1 = os.path.join(temp_dir, 'abc.tfa')
        fasta2 = os.path.join(temp_dir, 'def.tfa')
        
        try:
            with open(fasta1, 'w') as f:
                f.write('>abc_1\n')
                f.write('ATCGATCGATCGATCG\n')
            
            with open(fasta2, 'w') as f:
                f.write('>def_1\n')
                f.write('GCTAGCTAGCTAGCTA\n')
            
            # Load profile and fastas
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            fastas = Fastas(logger, temp_dir, 5, False)
            
            # Verify integration
            self.assertTrue(profile.has_gene('abc'))
            self.assertTrue(profile.has_gene('def'))
            self.assertGreater(len(fastas.fastas_to_kmers), 0)
            
        finally:
            # Cleanup
            for f in [profile_file, fasta1, fasta2]:
                if os.path.exists(f):
                    os.unlink(f)
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

    def test_fasta_and_kmers_integration(self):
        """Test that Fasta and Kmers modules work together"""
        logger = logging.getLogger(__name__)
        
        # Create a temporary fasta file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>seq1\n')
            f.write('ATCGATCGATCGATCGATCG\n')
            f.write('>seq2\n')
            f.write('GCTAGCTAGCTAGCTAGCTA\n')
            fasta_file = f.name
        
        try:
            fasta = Fasta(logger, fasta_file, 5, False)
            
            # Verify kmers were extracted
            self.assertIsInstance(fasta.all_kmers, dict)
            self.assertGreater(len(fasta.all_kmers), 0)
            
            # Verify sequences are present
            self.assertIn('seq1', fasta.sequences_to_kmers)
            self.assertIn('seq2', fasta.sequences_to_kmers)
            
        finally:
            if os.path.exists(fasta_file):
                os.unlink(fasta_file)

    def test_gene_coverage_calculation(self):
        """Test gene coverage calculation"""
        gene1 = Gene('ABC_1', 100, 0)
        gene2 = Gene('DEF_2', 50, 50)
        gene3 = Gene('GHI_3', 0, 100)
        
        self.assertTrue(gene1.is_full_coverage())
        self.assertFalse(gene2.is_full_coverage())
        self.assertFalse(gene3.is_full_coverage())
        
        self.assertEqual(str(gene1), 'ABC(1)')
        self.assertEqual(str(gene2), 'DEF(2)*')
        self.assertEqual(str(gene3), 'GHI(3)*')

    def test_kmer_extraction_from_sequence(self):
        """Test that k-mers are correctly extracted from sequences"""
        sequence = 'ATCGATCG'
        kmers = Kmers(sequence, 4)
        
        result = kmers.get_all_kmers_counter(max_kmer_count=10)
        
        # Verify result is a dictionary
        self.assertIsInstance(result, dict)
        
        # Expected k-mers: ATCG, TCGA, CGAT, GATC, ATCG (duplicate)
        expected_kmers = {'ATCG', 'TCGA', 'CGAT', 'GATC'}
        self.assertEqual(set(result.keys()), expected_kmers)

    def test_multiple_fasta_files_processing(self):
        """Test processing multiple FASTA files"""
        logger = logging.getLogger(__name__)
        
        # Create temporary directory with multiple fasta files
        temp_dir = tempfile.mkdtemp()
        fasta1 = os.path.join(temp_dir, 'file1.tfa')
        fasta2 = os.path.join(temp_dir, 'file2.tfa')
        
        try:
            with open(fasta1, 'w') as f:
                f.write('>gene1\n')
                f.write('ATCGATCGATCG\n')
            
            with open(fasta2, 'w') as f:
                f.write('>gene2\n')
                f.write('GCTAGCTAGCTA\n')
            
            fastas = Fastas(logger, temp_dir, 5, False)
            
            # Verify both files were processed
            self.assertEqual(len(fastas.filenames), 2)
            self.assertGreater(len(fastas.fastas_to_kmers), 0)
            
        finally:
            for f in [fasta1, fasta2]:
                if os.path.exists(f):
                    os.unlink(f)
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)

if __name__ == '__main__':
    unittest.main()
