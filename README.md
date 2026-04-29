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

Metadata files should follow qiime input conventions
- See : https://docs.qiime2.org/2024.10/tutorials/metadata/

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
```
Activate your QIIME2 environment:

```bash
conda activate qiime2-amplicon-2024.10
```
Make the script executable:

```bash
chmod +x qiimescraper.py
```
## FASTQ naming expectation

By default, QIIME Scraper expects paired-end FASTQ files to use Illumina-style names ending in:

```bash
_R1_001.fastq.gz
_R2_001.fastq.gz
```

Example:
```bash
Sample01_R1_001.fastq.gz
Sample01_R2_001.fastq.gz
Sample02_R1_001.fastq.gz
Sample02_R2_001.fastq.gz
```

The sample ID is inferred by removing the R1 suffix.

Custom suffixes can be supplied using:

```bash
--r1-suffix
--r2-suffix
```

## Primer input options

You can provide primers directly through the command line:

```bash
--forward-primers CTTGGTCATTTAGAGGAAGTAA \
--reverse-primers GCTGCGTTCTTCATCGATGC
```

For multiple primers, separate them with commas:
```bash
--forward-primers CTTGGTCATTTAGAGGAAGTAA,GCTGCGTTCTTCATCGATGC \
--reverse-primers GCTGCGTTCTTCATCGATGC,CTTGGTCATTTAGAGGAAGTAA
```

Alternatively, use a primer TSV file:

```bash
--primers-file primers.example.tsv --primer-set ITS
```

## Run Full Pipeline:

```bash
python qiimescraper.py \
  --full \
  --fastq-dir raw_fastqs \
  --outdir outputs/ITS_run \
  --primers-file primers.example.tsv \
  --primer-set ITS \
  --classifier classifiers/its_classifier.qza \
  --metadata metadata.tsv \
  --cores 8
```


### Run full pipeline with primers supplied directly:

```bash
python qiimescraper.py \
  --full \
  --fastq-dir raw_fastqs \
  --outdir outputs/my_run \
  --forward-primers CTTGGTCATTTAGAGGAAGTAA \
  --reverse-primers GCTGCGTTCTTCATCGATGC \
  --classifier classifiers/classifier.qza \
  --metadata metadata.tsv \
  --cores 8
```

### Main Outputs: 
```bash
outputs/my_run/
├── manifest.tsv
├── demux.qza
├── demux.qzv
├── demux_trimmed.qza
├── cutadapt_verbose.log
├── table.qza
├── repseqs.qza
├── dada2stats.qza
├── dada2stats.qzv
├── taxonomy.qza
├── taxonomy.qzv
└── taxa-barplot.qzv
```

## NOTES
- This tool does not replace QIIME2. It is a convenience wrapper around common QIIME2 commands.

- Always inspect demux.qzv before deciding DADA2 truncation parameters.

- Use a classifier appropriate for your marker region and database.

- For production or publication workflows, record your QIIME2 version, database version, primer sequences, and classifier source.

## Disclaimer
QIIME Scraper is intended as a helper utility for reproducible amplicon analysis. Users are responsible for validating parameters and outputs for their own datasets.update
update2
YOLO test
