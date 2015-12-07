#!/usr/bin/env python

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# Explicitly list bin scripts to be installed, seeing as I have a few local
# bin files that are not (yet) part of the distribution.
scripts = [
    'bin/adaptor-distances.py',
    'bin/alignments-per-read.py',
    'bin/cat-json-blast-records.py',
    'bin/check-fasta-json-blast-consistency.py',
    'bin/convert-blast-xml-to-json.py',
    'bin/convert-fasta-to-one-sequence-per-line.py',
    'bin/fasta-subset.py',
    'bin/fasta-subtraction.py',
    'bin/filter-fasta-by-length.py',
    'bin/filter-hits-to-fasta.py',
    'bin/find-hits.py',
    'bin/get-features.py',
    'bin/get-reads.py',
    'bin/graph-evalues.py',
    'bin/local-align.py',
    'bin/noninteractive-alignment-panel.py',
    'bin/pre-commit.sh',
    'bin/print-blast-xml-for-derek.py',
    'bin/print-blast-xml.py',
    'bin/print-read-lengths.py',
    'bin/read-blast-json.py',
    'bin/read-blast-xml.py',
    'bin/sff-to-fastq.py',
    'bin/split-fasta-by-adaptors.py',
    'bin/summarize-reads.py',
    'bin/trim-primers.py',
    'bin/trim-reads.py',
    'bin/write-htcondor-job-spec.py',
]

setup(name='dark-matter',
      version='1.0.1',
      packages=['dark'],
      url='https://github.com/acorg/dark-matter',
      download_url='https://github.com/acorg/dark-matter',
      author='Terry Jones, Barbara Muehlemann, Sophie Mathias',
      author_email='tcj25@cam.ac.uk',
      keywords=['virus discovery'],
      classifiers=[
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      license='MIT',
      description=('Python classes to help with virus discovery from genetic '
                   'sequence data'),
      scripts=scripts,
      ext_modules=cythonize([
          Extension('dark.gor4',
                    ['src/gor4/gor4.pyx', 'src/gor4/gor4-base.c',
                     'src/gor4/nrutil.c', 'src/gor4/api.c'])]
      ))
