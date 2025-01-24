# COVER

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Prepare Database](#prepare-database)
4. [Find Candidate Variant Pairs and Estimate Coverage](#find-candidate-variant-pairs-and-estimate-coverage)
   - [Input](#input)
   - [Output](#output)
   - [Parameters](#parameters)

## Overview

Code for COVER strategy. The document is under construction.

## Installation

Clone this repository `git clone https://github.com/Han-Cao/COVER.git`.

Install required tools:
- bcftools
- python3

Install python modules:

```
pip install numpy==1.26.4 pandas gtfparse pysam
```

## Prepare database

Run `setup.sh` to download ENSEMBL GTF annotation and 1000 Genomes VCF file and set up the transcirpt database. You will get:

- `data/Homo_sapiens.GRCh38.110.gtf.gz`: raw ENSEMBL GTF
- `data/Homo_sapiens.GRCh38.110.db`: SQLite database of the GTF
- `data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf`: High-coverage 1000 Genomes Project VCF file of unrelated individuals, only keep biallelic SNV with MAF > 0.01.

By default, ENSEMBL V110 is used. To use another version, modify `ENSEMBL_VERSION` in `setup.sh`.

## Find candidate variant pairs and estimate coverage

Run `run_all_het.py` to identify common variant pairs and estimate their coverage (i.e., heterozygous frequency) in a specific population.

```
python src/run_all_het.py \
-d data/Homo_sapiens.GRCh38.110.db \
-v data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf \
-l transcripts.txt \
-p EUR \
-o PREFIX
```

### Input

- GTF database: `data/Homo_sapiens.GRCh38.110.db`
- 1000 Genomes VCF: `data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf`
- Population: `AMR`, `AFR`, `EUR`, `EAS`, `SAS`
- Transcript list: one ENSEMBL transcript ID (without version suffix) per line. If not specified, all transcripts in the database will be processed.

### Output

- `PREFIX.het_freq.all.txt`: coverage of candidate variant pairs for each transcript
- `PREFIX.het_freq.top.txt`: variant pairs with the highest coverage for each transcript
- `PREFIX.pair_het_freq.all.txt`: coverage of the combination of 2 variant pairs for each transcript
- `PREFIX.pair_het_freq.top.txt`: combination of 2 variant pairs with the highest coverage for each transcript


### Parameters
```
run_all_het.py 

options:
  -d DB, --db DB        SQLite3 database file for GTF
  -v VCF, --vcf VCF     Reference 1000G VCF file
  -o OUTPUT, --output OUTPUT
                        Output file prefix
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
