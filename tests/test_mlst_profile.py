import unittest
import os
import tempfile
from krocus.MlstProfile import MlstProfile, Error

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data')

class TestMlstProfile(unittest.TestCase):

    def test_profile_initialization(self):
        """Test initialization with a valid profile file"""
        # Create a temporary profile file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\tghi\n')
            f.write('1\t1\t1\t1\n')
            f.write('2\t1\t1\t2\n')
            f.write('3\t2\t1\t1\n')
            profile_file = f.name
        
        try:
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            self.assertIsNotNone(profile)
            self.assertEqual(profile.genes_list, ['abc', 'def', 'ghi'])
            self.assertEqual(profile.genes_set, {'abc', 'def', 'ghi'})
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

    def test_has_gene(self):
        """Test has_gene method"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\tghi\n')
            f.write('1\t1\t1\t1\n')
            profile_file = f.name
        
        try:
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            self.assertTrue(profile.has_gene('abc'))
            self.assertTrue(profile.has_gene('def'))
            self.assertTrue(profile.has_gene('ghi'))
            self.assertFalse(profile.has_gene('xyz'))
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

    def test_get_sequence_type(self):
        """Test get_sequence_type method"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\tghi\n')
            f.write('1\t1\t1\t1\n')
            f.write('2\t1\t1\t2\n')
            f.write('3\t2\t1\t1\n')
            profile_file = f.name
        
        try:
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            # Test exact match
            self.assertEqual(profile.get_sequence_type({'abc': 1, 'def': 1, 'ghi': 1}), 1)
            self.assertEqual(profile.get_sequence_type({'abc': 1, 'def': 1, 'ghi': 2}), 2)
            self.assertEqual(profile.get_sequence_type({'abc': 2, 'def': 1, 'ghi': 1}), 3)
            # Test no match
            self.assertEqual(profile.get_sequence_type({'abc': 3, 'def': 3, 'ghi': 3}), 'ND')
            # Test incomplete match
            self.assertEqual(profile.get_sequence_type({'abc': 1, 'def': 1}), 'ND')
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

    def test_invalid_profile_file(self):
        """Test that error is raised for non-existent file"""
        with self.assertRaises(Error):
            MlstProfile('/nonexistent/file.txt')

    def test_invalid_header(self):
        """Test that error is raised for invalid header"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('WRONG\tabc\tdef\n')
            f.write('1\t1\t1\n')
            profile_file = f.name
        
        try:
            with self.assertRaises(Error):
                MlstProfile(profile_file)
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

    def test_columns_to_ignore(self):
        """Test that certain columns are ignored"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\tclonal_complex\tCC\tLineage\tmlst_clade\tspecies\n')
            f.write('1\t1\t1\t100\t200\tA\tB\tC\n')
            profile_file = f.name
        
        try:
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            # clonal_complex, CC, Lineage, mlst_clade, and species should be ignored
            self.assertEqual(profile.genes_list, ['abc', 'def'])
            self.assertNotIn('clonal_complex', profile.genes_list)
            self.assertNotIn('CC', profile.genes_list)
            self.assertNotIn('Lineage', profile.genes_list)
            self.assertNotIn('mlst_clade', profile.genes_list)
            self.assertNotIn('species', profile.genes_list)
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

    def test_duplicate_profiles_same_st(self):
        """Test handling of duplicate profiles with same ST"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\n')
            f.write('1\t1\t1\n')
            f.write('1\t1\t1\n')  # Duplicate
            profile_file = f.name
        
        try:
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            # Should handle duplicates gracefully
            self.assertEqual(profile.get_sequence_type({'abc': 1, 'def': 1}), 1)
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

    def test_duplicate_profiles_different_st(self):
        """Test handling of duplicate profiles with different STs"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('ST\tabc\tdef\n')
            f.write('1\t1\t1\n')
            f.write('2\t1\t1\n')  # Same alleles, different ST
            profile_file = f.name
        
        try:
            # With warnings disabled
            profile = MlstProfile(profile_file, duplicate_warnings=False)
            # Should use the smaller ST number
            result = profile.get_sequence_type({'abc': 1, 'def': 1})
            self.assertEqual(result, 1)
        finally:
            if os.path.exists(profile_file):
                os.unlink(profile_file)

if __name__ == '__main__':
    unittest.main()
