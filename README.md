# QIIME Scraper

QIIME Scraper is a lightweight command-line wrapper for running a paired-end QIIME2 amplicon workflow from raw FASTQ files to taxonomy barplots.

It was designed to make repetitive QIIME2 workflows easier to run while remaining flexible across marker genes such as 16S, ITS, COI, or custom amplicons.

## Features

- Generates QIIME2 paired-end manifest files automatically
- Imports paired-end FASTQ files into QIIME2
- Trims primer sequences with `qiime cutadapt trim-paired`
- Runs DADA2 paired-end denoising
- Assigns taxonomy using a trained QIIME2 classifier
- Creates interactive QIIME2 taxa barplots
- Supports custom primer sets through command-line arguments or a TSV primer file

## Requirements

This project assumes you already have a working QIIME2 environment installed.

Tested conceptually with a standard QIIME2 amplicon workflow using:

- QIIME2
- q2-cutadapt
- q2-dada2
- q2-feature-classifier
- Python 3.10+

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/qiimescraper.git
cd qiimescraper
