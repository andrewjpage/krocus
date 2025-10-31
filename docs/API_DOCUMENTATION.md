# Krocus API Documentation

## Overview

This document provides detailed API documentation for the Krocus Python modules. This is intended for developers who want to use Krocus as a library or contribute to the project.

## Core Modules

### Krocus

Main orchestration class that coordinates MLST analysis.

```python
from krocus.Krocus import Krocus
```

#### Class: `Krocus(options)`

**Parameters:**
- `options`: Namespace object with the following attributes:
  - `allele_directory` (str): Path to MLST database directory
  - `input_fastq` (str): Path to input FASTQ file
  - `kmer` (int): K-mer size (5-31)
  - `verbose` (bool): Enable debug logging
  - `min_fasta_hits` (int): Minimum k-mer matches
  - `print_interval` (int): Print results every N reads
  - `output_file` (str): Output file path (optional)
  - `filtered_reads_file` (str): Filtered reads output (optional)
  - `target_st` (int): Target ST for performance testing (optional)
  - `max_gap` (int): Maximum gap for blocks
  - `min_block_size` (int): Minimum block size in bases
  - `margin` (int): Flanking region size
  - `divisible_by_3` (bool): Filter genes not divisible by 3
  - `min_kmers_for_onex_pass` (int): Minimum k-mers for first pass
  - `max_kmers` (int): Maximum k-mer count threshold

**Methods:**
- `run()`: Execute MLST analysis
- `mlst_profile_file()`: Get path to profile.txt file

**Example:**
```python
from argparse import Namespace
from krocus.Krocus import Krocus

options = Namespace(
    allele_directory='db/',
    input_fastq='reads.fastq',
    kmer=11,
    verbose=False,
    min_fasta_hits=10,
    print_interval=500,
    output_file=None,
    filtered_reads_file=None,
    target_st=None,
    max_gap=4,
    min_block_size=150,
    margin=50,
    divisible_by_3=False,
    min_kmers_for_onex_pass=10,
    max_kmers=10
)

krocus = Krocus(options)
krocus.run()
```

---

### Kmers

K-mer extraction and counting.

```python
from krocus.Kmers import Kmers
```

#### Class: `Kmers(sequence, k)`

**Parameters:**
- `sequence` (str): DNA sequence
- `k` (int): K-mer size

**Methods:**
- `get_all_kmers_counter(max_kmer_count=1)`: Get k-mers as dict with 0 values
- `get_all_kmers_freq(max_kmer_count=10)`: Get k-mers with their frequencies
- `get_all_kmers(max_kmer_count=1)`: Get k-mers with KmerHit objects
- `get_one_x_coverage_of_kmers()`: Get k-mers with 1x coverage

**Returns:**
- Dictionary mapping k-mer strings to counts or KmerHit objects

**Example:**
```python
from krocus.Kmers import Kmers

sequence = "ATCGATCGATCG"
kmers = Kmers(sequence, k=5)

# Get k-mer frequencies
freq = kmers.get_all_kmers_freq(max_kmer_count=10)
# {'ATCGA': 2, 'TCGAT': 2, 'CGATC': 2, 'GATCG': 2}
```

---

### Fasta

Read and process FASTA files.

```python
from krocus.Fasta import Fasta
```

#### Class: `Fasta(logger, filename, k, divisible_by_3, max_kmers=5)`

**Parameters:**
- `logger`: Python logger object
- `filename` (str): Path to FASTA file
- `k` (int): K-mer size
- `divisible_by_3` (bool): Filter sequences not divisible by 3
- `max_kmers` (int): Maximum k-mer count threshold

**Attributes:**
- `sequences_to_kmers`: Dict mapping sequence IDs to k-mer counts
- `sequences_to_kmers_count`: Dict mapping sequence IDs to k-mer frequencies
- `all_kmers`: Dict of all k-mers across all sequences

**Methods:**
- `sequence_kmers()`: Extract k-mers from all sequences
- `sequence_kmers_vals()`: Extract k-mer frequencies
- `all_kmers_in_file()`: Count k-mers across file

**Example:**
```python
import logging
from krocus.Fasta import Fasta

logger = logging.getLogger(__name__)
fasta = Fasta(logger, 'gene.fa', k=11, divisible_by_3=False)

# Access k-mers
for seq_id, kmers in fasta.sequences_to_kmers.items():
    print(f"{seq_id}: {len(kmers)} k-mers")
```

---

### Fastas

Process multiple FASTA files.

```python
from krocus.Fastas import Fastas
```

#### Class: `Fastas(logger, allele_directory, k, divisible_by_3, max_kmers=10)`

**Parameters:**
- `logger`: Python logger object
- `allele_directory` (str): Directory containing .tfa files
- `k` (int): K-mer size
- `divisible_by_3` (bool): Filter sequences not divisible by 3
- `max_kmers` (int): Maximum k-mer count threshold

**Attributes:**
- `filenames`: List of FASTA file paths
- `fastas_to_kmers`: Dict mapping Fasta objects to k-mer dicts

**Methods:**
- `get_fastas_to_kmers()`: Process all FASTA files
- `allele_filenames(allele_directory)`: Find all .tfa files

**Example:**
```python
import logging
from krocus.Fastas import Fastas

logger = logging.getLogger(__name__)
fastas = Fastas(logger, 'allele_dir/', k=11, divisible_by_3=False)

print(f"Found {len(fastas.filenames)} allele files")
```

---

### Gene

Represent an allele with coverage information.

```python
from krocus.Gene import Gene
```

#### Class: `Gene(name, kmers_found, kmers_not_found)`

**Parameters:**
- `name` (str): Gene name (format: geneName_alleleNumber)
- `kmers_found` (int): Number of k-mers found
- `kmers_not_found` (int): Number of k-mers not found

**Methods:**
- `allele_name()`: Extract gene name
- `allele_number()`: Extract allele number
- `is_full_coverage()`: Check if all k-mers found
- `__str__()`: String representation (e.g., "geneX(1)" or "geneX(1)*")

**Example:**
```python
from krocus.Gene import Gene

gene = Gene('infB_1', kmers_found=100, kmers_not_found=0)
print(gene.allele_name())    # 'infB'
print(gene.allele_number())  # 1
print(gene.is_full_coverage())  # True
print(str(gene))  # 'infB(1)'

partial_gene = Gene('pgi_2', kmers_found=50, kmers_not_found=10)
print(str(partial_gene))  # 'pgi(2)*'
```

---

### Read

Represent a FASTQ read.

```python
from krocus.Read import Read
```

#### Class: `Read(id=None, seq=None, qual=None)`

**Parameters:**
- `id` (str): Read identifier
- `seq` (str): DNA sequence
- `qual` (str): Quality string

**Methods:**
- `subsequence(start, end)`: Extract subsequence
- `reverse_complement_sequence()`: Get reverse complement
- `reverse_read()`: Create reverse complement Read object
- `get_next_from_file(fh)`: Read next entry from FASTQ file
- `__str__()`: FASTQ format string

**Example:**
```python
from krocus.Read import Read

read = Read(id='read1', seq='ATCGATCG', qual='IIIIIIII')
print(str(read))

# Subsequence
subseq = read.subsequence(0, 4)
print(subseq.seq)  # 'ATCG'

# Reverse complement
rev_read = read.reverse_read()
print(rev_read.seq)  # 'CGATCGAT'
```

---

### Blocks

Find contiguous k-mer match blocks.

```python
from krocus.Blocks import Blocks
```

#### Class: `Blocks(k, max_gap, margin, min_block_size)`

**Parameters:**
- `k` (int): K-mer size
- `max_gap` (int): Maximum gap in k-mer multiples
- `margin` (int): Flanking region size
- `min_block_size` (int): Minimum block size in bases

**Methods:**
- `find_all_blocks(sequence_hits)`: Find all contiguous blocks
- `merge_blocks(blocks)`: Merge nearby blocks
- `find_largest_block(hits)`: Find largest merged block
- `adjust_block_start(start)`: Add margin to block start
- `adjust_block_end(end, seq_length)`: Add margin to block end

**Example:**
```python
from krocus.Blocks import Blocks

blocks = Blocks(k=11, max_gap=4, margin=50, min_block_size=150)
hits = [0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0]
blocks_found = blocks.find_all_blocks(hits)
# [[3, 14]]
```

---

### MlstProfile

Parse and query MLST profile.

```python
from krocus.MlstProfile import MlstProfile
```

#### Class: `MlstProfile(infile, duplicate_warnings=True)`

**Parameters:**
- `infile` (str): Path to profile.txt
- `duplicate_warnings` (bool): Print warnings for duplicate profiles

**Attributes:**
- `genes_list`: List of gene names (ordered)
- `genes_set`: Set of gene names
- `profile_to_type`: Dict mapping allele tuples to ST numbers

**Methods:**
- `has_gene(gene)`: Check if gene exists in profile
- `get_sequence_type(type_dict)`: Get ST from allele dictionary

**Example:**
```python
from krocus.MlstProfile import MlstProfile

profile = MlstProfile('profile.txt')

# Check gene
if profile.has_gene('infB'):
    print("infB is in the scheme")

# Get ST
alleles = {'infB': 1, 'pgi': 1, 'phoE': 9, 'tonB': 93}
st = profile.get_sequence_type(alleles)
print(f"ST: {st}")
```

---

### InputTypes

Input validation utilities.

```python
from krocus.InputTypes import InputTypes
```

**Static Methods:**
- `is_fastq_file_valid(filename)`: Validate FASTQ file exists
- `is_allele_directory_valid(filename)`: Validate directory exists
- `is_kmer_valid(value_str)`: Validate k-mer size (5-31)

**Example:**
```python
from krocus.InputTypes import InputTypes

# These are typically used with argparse
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fastq', type=InputTypes.is_fastq_file_valid)
parser.add_argument('--kmer', type=InputTypes.is_kmer_valid, default=11)
```

---

### KrocusDatabaseDownloader

Download MLST databases from PubMLST.

```python
from krocus.KrocusDatabaseDownloader import KrocusDatabaseDownloader
```

#### Class: `KrocusDatabaseDownloader(options)`

**Parameters:**
- `options`: Namespace with:
  - `list_species` (bool): List available species
  - `species` (str): Species name to download
  - `output_directory` (str): Output directory
  - `verbose` (bool): Enable verbose output

**Methods:**
- `run()`: Execute download operation

**Example:**
```python
from argparse import Namespace
from krocus.KrocusDatabaseDownloader import KrocusDatabaseDownloader

options = Namespace(
    list_species=False,
    species='Escherichia coli',
    output_directory='ecoli_db',
    verbose=True
)

downloader = KrocusDatabaseDownloader(options)
downloader.run()
```

---

## Testing

Krocus includes a comprehensive test suite with 78 tests and 65% code coverage.

### Running Tests

```bash
# Run all tests
python3 -m unittest discover -s tests/ -p 'test_*.py'

# Run with coverage
python3 -m coverage run -m unittest discover -s tests/
python3 -m coverage report --include='krocus/*.py'

# Run specific test
python3 -m unittest tests.test_kmers.TestKmers.test_four_kmers
```

### Test Organization

- `tests/test_blocks.py`: Blocks module tests
- `tests/test_gene.py`: Gene module tests
- `tests/test_kmers.py`: Kmers module tests
- `tests/test_read.py`: Read module tests
- `tests/test_fasta.py`: Fasta module tests
- `tests/test_fastas.py`: Fastas module tests
- `tests/test_fastq.py`: Fastq module tests
- `tests/test_krocus.py`: Main Krocus class tests
- `tests/test_mlst_profile.py`: MlstProfile tests
- `tests/test_input_types.py`: InputTypes tests
- `tests/test_krocus_database_downloader.py`: Database downloader tests
- `tests/test_integration.py`: Integration tests

### Writing New Tests

```python
import unittest
from krocus.YourModule import YourClass

class TestYourClass(unittest.TestCase):
    def test_something(self):
        obj = YourClass(params)
        result = obj.method()
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()
```

## Error Handling

All modules raise appropriate exceptions:

- `krocus.MlstProfile.Error`: Profile file errors
- `krocus.Fastq.Error`: FASTQ processing errors
- `krocus.PubmlstGetter.Error`: Download errors
- `argparse.ArgumentTypeError`: Input validation errors

## Logging

Krocus uses Python's logging module:

```python
import logging

# Enable debug logging
logger = logging.getLogger('krocus')
logger.setLevel(logging.DEBUG)
```

## Dependencies

- **biopython**: FASTA/FASTQ parsing
- **pyfastaq**: Sequence manipulation
- **numpy**: Numerical operations (Fastq module)

## Performance Considerations

1. **K-mer size**: Larger k-mers = faster but less sensitive
2. **max_kmers parameter**: Lower values = faster processing
3. **Memory usage**: ~1MB RAM per 1MB FASTQ
4. **I/O**: Streaming from stdin is more efficient than large files

## Contributing

See DEVELOPER_GUIDE.md for contribution guidelines.

## Version History

See CHANGELOG for version history and updates.
