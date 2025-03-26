#!/bin/bash

# specify ensemble version
ENSEMBL_VERSION=110

# download gtf from ensembl
wget \
-O data/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gtf.gz \
https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gtf.gz

# parse gtf and store into database
python3 src/gtf2db.py \
--gtf data/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gtf.gz \
--db data/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.db

# download 1000G phased vcf
# VCF requirement:
# fully phased, biallelic, no missing genotypes
[[ ! -d data/vcf_raw/ ]] && mkdir -p data/vcf_raw/
for i in {1..22} X; do
    wget \
    -O data/vcf_raw/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
done

# process vcf
[[ -e data/vcf_raw/vcf_list.txt ]] && rm data/vcf_raw/vcf_list.txt
for i in {1..22} X; do
    bcftools view \
    -m 2 -M 2 --type snps \
    -S data/1kg.unrelated.txt \
    -Ou data/vcf_raw/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
    bcftools annotate -x ^INFO/MAF_EUR_unrel,INFO/MAF_EAS_unrel,INFO/MAF_AMR_unrel,INFO/MAF_SAS_unrel,INFO/MAF_AFR_unrel,INFO/ExcHet_EAS,INFO/ExcHet_EUR,INFO/ExcHet_AMR,INFO/ExcHet_SAS,INFO/ExcHet_AFR \
    -Ob -o data/vcf_raw/1kGP_high_coverage_Illumina.chr${i}.filtered.unrelated.SNV.MAF1.bcf
    echo "data/vcf_raw/1kGP_high_coverage_Illumina.chr${i}.filtered.unrelated.SNV.MAF1.bcf" >> data/vcf_raw/vcf_list.txt
done

# concat vcf
bcftools concat -n -f data/vcf_raw/vcf_list.txt -o data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf
bcftools index data/1kGP_high_coverage_Illumina.filtered.unrelated.SNV.MAF1.bcf
