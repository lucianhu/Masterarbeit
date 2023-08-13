#!/bin/bash

# Go to SnpEff's install dir
cd $HOME/DNA_softwares/snpEff

# Create database dir
mkdir data/CHM13.2
cd data/CHM13.2

# Download UCSC human reference genome
aria2c -x8 https://hgdownload-test.gi.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz -o sequences.fa.gz
gzip -d sequences.fa.gz
samtools faidx sequences.fa

# Download Ensemble annotated genes
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz -o genes.gtf.gz

# Download Ensemble proteins
# This is used for:
#   - "Rare Amino Acid" annotations
#   - Sanity check (checking protein predicted from DNA sequences match 'real' proteins)
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-pep.fa.gz -o protein.fa.gz

# Download Ensemble transcripts
# Note: This is used as "sanity check" (checking that CDSs predicted from gene sequences match 'real' CDSs)
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-cdna.fa.gz -o cds.fa.gz

# Uncompress
gunzip *.gz

# Download Ensemble-ClinVar

aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-clinvar.vcf.gz

## Fix ClinVar file
zcat Homo_sapiens-GCA_009914755.4-2022_10-clinvar.vcf.gz | sed 's/##INFO=<ID=T2T-CHM13v2.0.gff.gz,Number=.,Type=String,Description=""/##INFO=<ID=T2T-CHM13v2.0.gff.gz,Number=.,Type=String,Description="UCSC HPRC Assembly Hub">/' | bgzip -c > CHM13.2.clinvar.vcf.gz

gzip -d CHM13.2.clinvar.vcf.gz

bgzip CHM13.2.clinvar.vcf && tabix -p vcf CHM13.2.clinvar.vcf.gz

# Uncompress:
# Why do we need to uncompress?
# Because ENSEMBL compresses files using a block compress gzip which is not compatible with Java's library Gunzip

# Edit snpEff.config file
#
# WARNING! You must do this yourself. Just copying and pasting this into a terminal won't work.
#
# Add lines:
# T2T-CHM13v2.0 release from ENSEMBL
# CHM13.2.genome : Human genome (Homo sapiens) T2T-CHM13v2.0 using transcripts
# CHM13.2.reference: https://projects.ensembl.org/hprc/

# Now we are ready to build the database
java -Xmx20g -jar $HOME/DNA_softwares/snpEff/snpEff.jar build -v CHM13.2 2>&1 | tee CHM13.2.build
