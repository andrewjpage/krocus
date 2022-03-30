import os
import shutil
import sys
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = 'x.y.z'
if os.path.exists('VERSION'):
  version = open('VERSION').read().strip()

setup(
    name='krocus',
    version=version,
    description='krocus: multi-locus sequence typing from uncorrected long reads',
	long_description='krocus: multi-locus sequence typing from uncorrected long reads',
    packages = find_packages(),
    author='Andrew J. Page',
    author_email='andrewjpage+krocus@gmail.com',
    url='https://github.com/andrewjpage/krocus',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
           'biopython >= 1.68',
		   'pyfastaq >= 3.12.0'
       ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)