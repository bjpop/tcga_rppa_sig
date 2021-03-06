#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''The program reads one or more input FASTA files.
For each file it computes a variety of statistics, and then
prints a summary of the statistics as output.

The goal is to provide a solid foundation for new bioinformatics command line tools,
and is an ideal starting place for new projects.'''


setup(
    name='tcga_rppa_sig',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['tcga_rppa_sig'],
    package_dir={'tcga_rppa_sig': 'tcga_rppa_sig'},
    entry_points={
        'console_scripts': ['tcga_rppa_sig = tcga_rppa_sig.tcga_rppa_sig:main']
    },
    url='https://github.com/bjpop/tcga_rppa_sig',
    license='LICENSE',
    description=('A prototypical bioinformatics command line tool'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython"],
)
