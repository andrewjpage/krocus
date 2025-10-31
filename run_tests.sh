#!/bin/bash
# Run Krocus test suite
# Tests are located in the tests/ directory

python3 -m unittest discover -s tests/ -p 'test_*.py'

