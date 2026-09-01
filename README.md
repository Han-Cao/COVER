# COVER

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Prepare Database](#prepare-database)
4. [Find Candidate Variant Pairs and Estimate Coverage](#find-candidate-variant-pairs-and-estimate-coverage)
   - [Input](#input)
   - [Output](#output)
   - [Parameters](#parameters)
5. [COVER API](#cover-api)

## Overview

Code for COVER strategy. The document is under construction.

## Installation

COVER was developed in Python >= 3.10.

To install COVER, clone this repository `git clone https://github.com/Han-Cao/COVER.git` and run:

```
pip install .
```

To use COVER-api, install additional dependencies:

```
pip install .[api]
```

To prepare required 1000 Genomes VCFs, make sure you have `bcftools` installed.

## Prepare database

Run `setup.sh` to download 1000 Genomes VCF file. You will get:

- `data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf`: High-coverage 1000 Genomes Project VCF file of unrelated individuals, only keep biallelic SNV with MAF > 0.01.

Transcript annotation of ENSEMBL V110 is provided in `data/Homo_sapiens.GRCh38.110.db`. To use another ENSEMBL version, run `cover-gtf2db`:

```
usage: cover-gtf2db [-h] [-g GTF] [-d DB]

Parse GTF files and store into a database

options:
  -h, --help         show this help message and exit
  -g GTF, --gtf GTF  input GTF file
  -d DB, --db DB     output SQLite3 database file
```

## Find candidate variant pairs and estimate coverage

Run `cover` to identify common variant pairs and estimate their coverage (i.e., heterozygous frequency) in a specific population.

```
cover \
-d data/Homo_sapiens.GRCh38.110.db \
-v data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf \
-l transcripts.txt \
-p EUR \
-o PREFIX \
--cpu 4
```

### Input

- GTF database: `data/Homo_sapiens.GRCh38.110.db`
- 1000 Genomes VCF: `data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf`
- Population: `AMR`, `AFR`, `EUR`, `EAS`, `SAS`
- Transcript list (optional): one ENSEMBL transcript ID (without version suffix) per line. If not specified, all transcripts in the database will be processed.
- Number of CPUs (optional): `4` (this is also the default value)

### Output

- `PREFIX.het_freq.all.txt`: coverage of candidate variant pairs for each transcript
- `PREFIX.het_freq.top.txt`: variant pairs with the highest coverage for each transcript
- `PREFIX.pair_het_freq.all.txt`: coverage of the combination of 2 variant pairs for each transcript
- `PREFIX.pair_het_freq.top.txt`: combination of 2 variant pairs with the highest coverage for each transcript


### Parameters
```
usage: cover [-h] -d DB -v VCF -o OUTPUT [--target-exons TARGET_EXONS] -p {AFR,AMR,EAS,EUR,SAS} [--cpu CPU] [-l ID_LIST] [--save-pickle] [--max-deletion MAX_DELETION] [--splice-donor-len SPLICE_DONOR_LEN]
             [--splice-receptor-len SPLICE_RECEPTOR_LEN] [--n-before-stop N_BEFORE_STOP] [--include-start-loss] [--start-loss-only] [--maf MAF] [--exchet EXCHET] [--pair-per-tx PAIR_PER_TX]
             [--n-pair-max N_PAIR_MAX] [--pair-het-cutoff PAIR_HET_CUTOFF] [--top-n-comb TOP_N_COMB]

Calculate co-heterozygous frequency for all transcripts in the database

options:
  -h, --help            show this help message and exit
  -d DB, --db DB        SQLite3 database file for GTF
  -v VCF, --vcf VCF     Reference 1000G VCF file
  -o OUTPUT, --output OUTPUT
                        Output file prefix
  --target-exons TARGET_EXONS
                        Predefined exons to target (default: None)
  -p {AFR,AMR,EAS,EUR,SAS}, --pop {AFR,AMR,EAS,EUR,SAS}
                        Population
  --cpu CPU             Number of CPUs to use
  -l ID_LIST, --id-list ID_LIST
                        List of transcript IDs to process
  --save-pickle         Save transcript db to pickle file
  --max-deletion MAX_DELETION
                        Maximum deletion size (default: 10000)
  --splice-donor-len SPLICE_DONOR_LEN
                        Length of splice donor region (default: 10)
  --splice-receptor-len SPLICE_RECEPTOR_LEN
                        Length of splice receptor region (default: 28)
  --n-before-stop N_BEFORE_STOP
                        Minimum number of exons before the stop codon to be considered as target (default: 2)
  --include-start-loss  Include start loss when finding target region
  --start-loss-only     Only include start loss when finding target region
  --maf MAF             MAF cutoff (default: 0.05)
  --exchet EXCHET       Excess heterozygosity test p-value cutoff (default: 1e-5)
  --pair-per-tx PAIR_PER_TX
                        Maximum number of variant pairs per transcript to output (default: 1000)
  --n-pair-max N_PAIR_MAX
                        Maximum number of variant pairs to be considered for combinations of two variant pairs (default: 200)
  --pair-het-cutoff PAIR_HET_CUTOFF
                        Minimum heterozygous frequency to be considered for combinations of two variant pairs (default: 0.1)
  --top-n-comb TOP_N_COMB
                        Top N combination of two variant pairs to be output (default: 10)
```

## COVER API

A web service of COVER is under development. The API server is provided in `cover-api`.

Start the API server:
```
cover-api \
-t data/Homo_sapiens.GRCh38.110.db \
-r data/results.db \
-v data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf \
--port 8000 \
--host 127.0.0.1
```

The `data/results.db` file is generated from the result file `PREFIX.het_freq.all.txt`:
```
cover-result2db -d results.db --het-freq-file PREFIX.het_freq.all.txt
```