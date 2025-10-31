# Krocus Developer Guide

## Overview

This guide is for developers who want to contribute to Krocus or understand its internal architecture.

## Architecture

Krocus follows a modular architecture with clear separation of concerns:

```
┌─────────────────────────────────────────┐
│          scripts/krocus                 │  (CLI Entry Point)
└─────────────────┬───────────────────────┘
                  │
         ┌────────▼─────────┐
         │  Krocus.py       │  (Main Orchestrator)
         └────────┬─────────┘
                  │
    ┌─────────────┼─────────────┐
    │             │             │
┌───▼───┐    ┌───▼───┐    ┌───▼────┐
│Fastas │    │Fastq  │    │ Mlst   │
│       │    │       │    │Profile │
└───┬───┘    └───┬───┘    └────────┘
    │            │
┌───▼───┐    ┌───▼───┐
│Fasta  │    │Read   │
│       │    │       │
└───┬───┘    └───┬───┘
    │            │
┌───▼────────────▼───┐
│      Kmers         │
└────────────────────┘
```

## Module Responsibilities

### Entry Points (scripts/)

- **krocus**: Main CLI for MLST analysis
- **krocus_database_downloader**: Download databases from PubMLST

### Core Analysis (krocus/)

- **Krocus.py**: Main orchestration and workflow management
- **Fastq.py**: Read processing and mapping
- **Fastas.py**: Multiple FASTA file management
- **Fasta.py**: Single FASTA file processing
- **MlstProfile.py**: MLST profile parsing and ST determination
- **Read.py**: FASTQ read representation
- **Kmers.py**: K-mer extraction and counting
- **Gene.py**: Allele representation with coverage
- **Blocks.py**: Contiguous region detection
- **InputTypes.py**: Input validation
- **KrocusDatabaseDownloader.py**: Database download orchestration
- **PubmlstGetter.py**: PubMLST API interaction

## Development Setup

### Clone Repository

```bash
git clone https://github.com/andrewjpage/krocus.git
cd krocus
```

### Create Virtual Environment

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### Install in Development Mode

```bash
pip install -e .
```

### Install Development Dependencies

```bash
pip install coverage
```

## Code Style

### Python Style

Follow PEP 8 guidelines:

- 4 spaces for indentation (no tabs)
- Maximum line length: 120 characters (flexible)
- Use descriptive variable names
- Add docstrings to classes and public methods

### Naming Conventions

- **Classes**: PascalCase (e.g., `MlstProfile`)
- **Functions/Methods**: snake_case (e.g., `get_sequence_type`)
- **Constants**: UPPER_SNAKE_CASE (e.g., `MAX_KMER_SIZE`)
- **Private methods**: Leading underscore (e.g., `_internal_method`)

### Documentation

- Use docstrings for all public classes and methods
- Include parameter types and return types
- Provide usage examples for complex functions

Example:
```python
def get_sequence_type(self, type_dict):
    """
    Determine sequence type from allele dictionary.
    
    Args:
        type_dict (dict): Dictionary mapping gene names to allele numbers
        
    Returns:
        int or str: Sequence type number, or 'ND' if not determined
        
    Example:
        >>> profile = MlstProfile('profile.txt')
        >>> st = profile.get_sequence_type({'gene1': 1, 'gene2': 5})
        >>> print(st)
        42
    """
```

## Testing

### Test Structure

Tests are located in `tests/` directory:

```
tests/
├── __init__.py
├── test_blocks.py
├── test_fasta.py
├── test_fastas.py
├── test_fastq.py
├── test_gene.py
├── test_input_types.py
├── test_integration.py
├── test_kmers.py
├── test_krocus.py
├── test_krocus_database_downloader.py
├── test_mlst_profile.py
├── test_read.py
└── data/
    ├── fasta/
    ├── fastas/
    └── fastq/
```

### Writing Tests

Use Python's unittest framework:

```python
import unittest
from krocus.YourModule import YourClass

class TestYourClass(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.obj = YourClass(params)
    
    def tearDown(self):
        """Clean up after tests"""
        pass
    
    def test_basic_functionality(self):
        """Test basic functionality"""
        result = self.obj.method()
        self.assertEqual(result, expected)
    
    def test_edge_case(self):
        """Test edge cases"""
        with self.assertRaises(ValueError):
            self.obj.method(invalid_input)
```

### Running Tests

```bash
# Run all tests
python3 -m unittest discover -s tests/ -p 'test_*.py'

# Run specific test file
python3 -m unittest tests.test_kmers

# Run specific test class
python3 -m unittest tests.test_kmers.TestKmers

# Run specific test method
python3 -m unittest tests.test_kmers.TestKmers.test_four_kmers

# Run with verbose output
python3 -m unittest discover -s tests/ -p 'test_*.py' -v
```

### Coverage Analysis

```bash
# Run tests with coverage
python3 -m coverage run -m unittest discover -s tests/

# Generate report
python3 -m coverage report --include='krocus/*.py'

# Generate HTML report
python3 -m coverage html --include='krocus/*.py'
# Open htmlcov/index.html in browser

# Check specific module
python3 -m coverage report --include='krocus/Kmers.py' -m
```

### Test Guidelines

1. **Test one thing per test method**
2. **Use descriptive test names**: `test_kmers_with_empty_sequence`
3. **Test edge cases**: empty inputs, boundary values, invalid inputs
4. **Test error conditions**: Use `assertRaises`
5. **Use test fixtures**: Create reusable test data
6. **Mock external dependencies**: Network calls, file I/O
7. **Keep tests fast**: Tests should run in seconds, not minutes
8. **Tests should be independent**: No shared state between tests

## Adding New Features

### 1. Plan the Feature

- Document the requirements
- Design the API
- Consider backward compatibility
- Plan test coverage

### 2. Implement the Feature

```python
# krocus/NewModule.py

class NewModule:
    """
    Brief description of the module.
    
    Detailed description of what it does and how it works.
    """
    
    def __init__(self, param1, param2):
        """
        Initialize the module.
        
        Args:
            param1 (type): Description
            param2 (type): Description
        """
        self.param1 = param1
        self.param2 = param2
    
    def new_method(self):
        """
        Method description.
        
        Returns:
            type: Description of return value
        """
        # Implementation
        pass
```

### 3. Write Tests

```python
# tests/test_new_module.py

import unittest
from krocus.NewModule import NewModule

class TestNewModule(unittest.TestCase):
    def test_initialization(self):
        """Test module initialization"""
        module = NewModule(param1='value1', param2='value2')
        self.assertEqual(module.param1, 'value1')
    
    def test_new_method(self):
        """Test new method functionality"""
        module = NewModule(param1='value1', param2='value2')
        result = module.new_method()
        self.assertIsNotNone(result)
```

### 4. Update Documentation

- Add API documentation to `docs/API_DOCUMENTATION.md`
- Add usage examples to `docs/USER_GUIDE.md`
- Update README.md if necessary

### 5. Run Tests and Coverage

```bash
python3 -m unittest discover -s tests/
python3 -m coverage run -m unittest discover -s tests/
python3 -m coverage report --include='krocus/*.py'
```

## Code Review Checklist

Before submitting a pull request:

- [ ] All tests pass
- [ ] New code has test coverage (aim for >80%)
- [ ] Code follows PEP 8 style guidelines
- [ ] Docstrings added for new functions/classes
- [ ] Documentation updated
- [ ] No unnecessary dependencies added
- [ ] Backward compatibility maintained
- [ ] Error handling implemented
- [ ] Edge cases considered

## Debugging

### Enable Verbose Logging

```bash
# CLI
krocus database_dir reads.fastq -v

# In code
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Use Python Debugger

```python
import pdb

def problematic_function():
    # ... some code ...
    pdb.set_trace()  # Breakpoint here
    # ... more code ...
```

### Profile Performance

```python
import cProfile
import pstats

cProfile.run('krocus.run()', 'profile_output')
p = pstats.Stats('profile_output')
p.sort_stats('cumulative').print_stats(20)
```

## Common Development Tasks

### Adding a New Command-Line Option

1. Edit `scripts/krocus`:
```python
parser.add_argument('--new_option', help='Description', type=int, default=10)
```

2. Update `Krocus.__init__()`:
```python
self.new_option = options.new_option
```

3. Use the option in your code

4. Update documentation

### Fixing a Bug

1. Write a failing test that demonstrates the bug
2. Fix the bug
3. Verify the test passes
4. Check that existing tests still pass
5. Update documentation if behavior changed

### Improving Performance

1. Profile the code to identify bottlenecks
2. Optimize the slow parts
3. Verify correctness with tests
4. Measure improvement
5. Document any changes

## Continuous Integration

Krocus uses Travis CI for automated testing:

- Tests run automatically on all commits and pull requests
- Coverage reports generated
- Build status visible on README

Configuration in `.travis.yml`

## Release Process

1. Update `VERSION` file
2. Update `CHANGELOG`
3. Run full test suite
4. Create git tag: `git tag -a v1.0.0 -m "Release 1.0.0"`
5. Push tag: `git push origin v1.0.0`
6. Build and upload to PyPI:
   ```bash
   python3 setup.py sdist
   twine upload dist/*
   ```

## Project Structure

```
krocus/
├── krocus/              # Main Python package
│   ├── __init__.py
│   ├── Blocks.py
│   ├── Fasta.py
│   ├── Fastas.py
│   ├── Fastq.py
│   ├── Gene.py
│   ├── InputTypes.py
│   ├── Kmers.py
│   ├── Krocus.py
│   ├── KrocusDatabaseDownloader.py
│   ├── MlstProfile.py
│   ├── PubmlstGetter.py
│   └── Read.py
├── scripts/             # Command-line scripts
│   ├── krocus
│   └── krocus_database_downloader
├── tests/               # Test suite
│   ├── test_*.py
│   └── data/
├── docs/                # Documentation
│   ├── USER_GUIDE.md
│   ├── API_DOCUMENTATION.md
│   ├── DEVELOPER_GUIDE.md
│   └── WORK_INSTRUCTION.md
├── databases/           # Bundled MLST databases
├── setup.py             # Package setup
├── README.md            # Project overview
├── LICENSE              # GPL v3 license
├── VERSION              # Version number
├── CHANGELOG            # Version history
└── .travis.yml          # CI configuration
```

## Resources

- **GitHub Repository**: https://github.com/andrewjpage/krocus
- **PubMLST**: https://pubmlst.org/
- **PeerJ Paper**: https://doi.org/10.7717/peerj.5233
- **Python Documentation**: https://docs.python.org/3/
- **BioPython**: https://biopython.org/

## Getting Help

- **Issues**: Report bugs or request features on GitHub Issues
- **Discussions**: Use GitHub Discussions for questions
- **Email**: andrewjpage+krocus@gmail.com

## License

Krocus is licensed under GNU GPL version 3. See LICENSE file for details.

When contributing, you agree that your contributions will be licensed under the same license.
