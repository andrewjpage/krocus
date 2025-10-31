import unittest
import os
import tempfile
from krocus.Read import Read

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data')

class TestRead(unittest.TestCase):

    def test_initialise(self):
        r = Read(id='read1', seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT', qual='I'*40)
        self.assertEqual(len(r.seq), 40)

    def test_one_read(self):
        r = Read(id='read1', seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT', qual='I'*40)
        subseq = r.subsequence(0, 4)
        self.assertEqual(subseq.seq, 'ACGT')
        subseq = r.subsequence(4, 8)
        self.assertEqual(subseq.seq, 'ACGT')

    def test_subsequence(self):
        r = Read(id='read1', seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT', qual='I'*40)
        self.assertEqual(r.subsequence(10, 20).seq, 'GTACGTACGT')
        self.assertEqual(r.subsequence(0, 40).seq, 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT')

    def test_subsequence_boundaries(self):
        r = Read(id='read1', seq='ACGTACGT', qual='I'*8)
        self.assertEqual(r.subsequence(0, 4).seq, 'ACGT')
        self.assertEqual(r.subsequence(4, 8).seq, 'ACGT')

    def test_read_name(self):
        r = Read(id='test_read_name', seq='ACGTACGT', qual='I'*8)
        self.assertEqual(r.id, 'test_read_name')

    def test_reverse_complement(self):
        r = Read(id='read1', seq='ATCG', qual='IIII')
        rev_comp = r.reverse_complement_sequence()
        self.assertEqual(rev_comp, 'CGAT')

    def test_reverse_read(self):
        r = Read(id='read1', seq='ATCG', qual='IIII')
        rev_read = r.reverse_read()
        self.assertEqual(rev_read.seq, 'CGAT')
        self.assertEqual(rev_read.id, 'read1_reverse')

    def test_str_representation(self):
        r = Read(id='read1', seq='ACGT', qual='IIII')
        output = str(r)
        self.assertIn('@read1', output)
        self.assertIn('ACGT', output)
        self.assertIn('+', output)
        self.assertIn('IIII', output)

    def test_get_next_from_file(self):
        """Test reading from a FASTQ file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write('@read1\n')
            f.write('ACGTACGT\n')
            f.write('+\n')
            f.write('IIIIIIII\n')
            f.write('@read2\n')
            f.write('GCTAGCTA\n')
            f.write('+\n')
            f.write('HHHHHHHH\n')
            fastq_file = f.name
        
        try:
            with open(fastq_file, 'r') as fh:
                r = Read()
                # Read first record
                result = r.get_next_from_file(fh)
                self.assertTrue(result)
                self.assertEqual(r.id, 'read1')
                self.assertEqual(r.seq, 'ACGTACGT')
                self.assertEqual(r.qual, 'IIIIIIII')
                
                # Read second record
                result = r.get_next_from_file(fh)
                self.assertTrue(result)
                self.assertEqual(r.id, 'read2')
                self.assertEqual(r.seq, 'GCTAGCTA')
                self.assertEqual(r.qual, 'HHHHHHHH')
                
                # No more records
                result = r.get_next_from_file(fh)
                self.assertFalse(result)
        finally:
            if os.path.exists(fastq_file):
                os.unlink(fastq_file)

    def test_get_next_from_file_with_blank_lines(self):
        """Test reading from a FASTQ file with blank lines"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write('\n')
            f.write('@read1\n')
            f.write('ACGTACGT\n')
            f.write('+\n')
            f.write('IIIIIIII\n')
            fastq_file = f.name
        
        try:
            with open(fastq_file, 'r') as fh:
                r = Read()
                result = r.get_next_from_file(fh)
                self.assertTrue(result)
                self.assertEqual(r.id, 'read1')
                self.assertEqual(r.seq, 'ACGTACGT')
        finally:
            if os.path.exists(fastq_file):
                os.unlink(fastq_file)

if __name__ == '__main__':
    unittest.main()

if __name__ == '__main__':
    unittest.main()
