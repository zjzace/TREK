#!/usr/bin/env python3
"""
Setup script for TREK
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    long_description = readme_file.read_text()
else:
    long_description = "TREK - Transcription End site identifier from long-read data"

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    requirements = requirements_file.read_text().strip().split('\n')
else:
    requirements = [
        'numpy>=1.19.0',
        'pandas>=1.2.0',
        'scikit-learn>=0.24.0',
        'pysam>=0.16.0',
        'biopython>=1.79',
        'interlap>=0.2.7',
        'joblib>=1.0.0',
        'tqdm>=4.60.0'
    ]

setup(
    name='TREK',
    version='1.0.0',
    description='TREK - Transcription End site identifier from long-read data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Jizhou ZHANG',
    author_email='jizhouzhang@cuhk.edu.hk',
    url='https://github.com/zjzace/TREK',
    packages=find_packages(),
    install_requires=requirements,
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'trek=trek.run_trek:main',
        ],
    },
    package_data={
        'trek': ['*.pl', '*.txt', '*.md'],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='bioinformatics RNA-seq polyadenylation long-read nanopore pacbio',
    project_urls={
        'Bug Reports': 'https://github.com/yourusername/polyAFinder/issues',
        'Source': 'https://github.com/yourusername/polyAFinder',
    },
)
