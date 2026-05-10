# CRISPRDataCraft

CRISPRDataCraft is a software pipeline designed to enable the identification and statistical analysis of genome editing events from next-generation sequencing data.

## Contents

    - Overview
    - Repo Contents
    - System Requirements
      - Hardware Requirements
      - Software Requirements
    - Installation Guide
      - Development Version
    - Usage
      - CRISPRDataCraft Main Pipeline
      - Codon Shift Analysis
      - Desired Mutant Analysis
    - Citation

## Overview

CRISPRDataCraft provides a complete workflow for analyzing genome editing outcomes from paired-end NGS data. The pipeline follows five key steps:

1. Barcode-based Read Splitting: NGS reads are aplit based on barcode sequences provide in ./example/barcode.csv
2. Paired-end Read Merging: Paired-end reads are merged using FLASH2 software (available at https://ccb.jhu.edu/software/FLASH/)
3. Target Sequence Extraction: The target sequence, flanked by two predefined fragments, is extracted using a reference sequence and flanking fragments provided in ./example/gene.fa
4. Customized Sequence Alignment: The target segments in the sequencing reads are aligned to the reference sequence using an customized algorithm. This alignment algorithm employs a "divide and conquer" strategy: it identifies matching seeds, extends them bidirectionally until mismatches are enountered, and recursively searches for new seeds in the remaining regions, ensuring both efficiency and accuracy in the alignment process.
5. Mutation Quantification: Mutations, insertions, and deletions are quantified to determine whether a read is modified or unmodified. This module was developed by modifying CRISPResso2 software (https://github.com/pinellolab/CRISPResso2/tree/master).
  
In addition, CRISPRDataCraft provides two downstream analysis tools using the output of 'IGDB_split_merge_seed_aligner_local.py' as input:
    - Codon Shift Analysis: Statistical evaluation of codon shifts resulting from genome editing.
      - Note: The reference fragment in 'gene.fa' must be a CDS starting at the first position of the protein-encoding triplet codon. 
    - Desired Mutant Analysis: Statistical assessment of desired mutants generated through genome editing.

---

## Repo Contents

- 'IGDB_split_merge_seed_aligner_local.py' – Main pipeline script.
- 'IGDB_analyze_codon_shift.py' – Codon shift analysis tool.
- 'IGDB_take_desired_mutants.py' – Desired mutant analysis tool.
- 'setup_compare.py' – Cython compilation script.
- 'example/' – Example input files (barcode.csv, gene.fa, clean fq.gz files).
- 'output_exp1/', 'output_exp2/' – Example output directories.

---
## System Requirements

### Hardware Requirements

The ‘CRISPRDataCraft’ package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 4 GB of RAM. For optimal performance, we recommend a computer with the following specs:

- RAM: 8+ GB
- CPU: 4+ cores, 2.5+ GHz/core

### Software Requirements

#### OS Requirements

The package development version is tested on Linux operating systems. The developmental version of the package has been tested on the following systems:

- Linux: Ubuntu 22.04

Before setting up the ‘CRISPRDataCraft’ package, users should have Python 3.12.7 or higher, and several packages installed from PyPI.

##### Installing Python 3.12.7 on Ubuntu 22.04

The latest version of Python can be installed by adding the deadsnakes repository to 'apt':

sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.12 python3.12-venv python3.12-dev

which should install in about 30 seconds.

## Installation Guide

### Development Version

#### Package Dependencies

Users should install the following packages prior to installing CRISPRDataCraft, from a terminal:
pip install numpy pandas biopython matplotlib seaborn plotly cython setuptools
which will install in about 1 minute on a machine with the recommended specs.

The 'CRISPRDataCraft' package functions with all packages in their latest versions as they appear on PyPI. The versions of software tested are, specifically:
Package	Version
numpy	≥ 1.24
pandas	≥ 2.0
biopython ≥ 1.81
matplotlib  ≥ 3.7
seaborn	    ≥ 0.12
plotly	    ≥ 5.14
cython	    ≥ 3.0
setuptools  ≥ 65.0
FLASH	    v2.2.00 (external binary)

If you are having an issue that you believe to be tied to software versioning, please open an Issue.
#### Package Compilation
python setup_compare.py build_ext --inplace

## Usage

### CRISPRDataCraft Main Pipeline
  
python3 IGDB_split_merge_seed_aligner_local.py -r1 ./example/1_clean_r1.fq.gz -r2 ./example/1_clean_r2.fq.gz -b ./example/barcode_1.csv -a ./example/gene_1.fa -p 5 -o output_exp1
  
python3 IGDB_split_merge_seed_aligner_local.py -r1 ./example/2_clean_r1.fq.gz -r2 ./example/2_clean_r2.fq.gz -b ./example/barcode_2.csv -a ./example/gene_2.fa -p 5 -o output_exp2


### Codon Shift Analysis:
  
python IGDB_analyze_codon_shift.py -pa output_exp1 -o1 ./output_exp1/CSAout1.csv -o2 ./output_exp1/CSAout2.csv  


### Desired Mutant Analysis:
  
python IGDB_take_desired_mutants.py -aft ./output_exp2/ -dsf ./example/ListbyProduct.txt -o ./output_exp2/Desired_mutants.csv -top 10 -pct 1
