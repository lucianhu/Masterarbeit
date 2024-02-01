#!/usr/bin/env bash
# Processing the human reference genome

# Download the human reference genome file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Unzip the downloaded file and save it as GRCh38.fa
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38.fa

# Download the exome file
wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED -O GRCh38-exome-interval.bed 

# Index the reference genome with samtools
samtools faidx GRCh38.fa

# Extract the introns/non-exome intervals using bedtools complement
bedtools complement -i GRCh38-exome-interval.bed -g GRCh38.fa.fai > GRCh38-intron-interval.bed

# Create a masked genome file by replacing intron/non-exome regions with "N" and retain only exons
bedtools maskfasta -fi GRCh38.fa -bed GRCh38-intron-interval.bed -fo GRCh38.exome.fa

# Index the masked reference genome with samtools
samtools faidx GRCh38.exome.fa

# Create a sequence dictionary for the reference using GATK
gatk CreateSequenceDictionary --R GRCh38.exome.fa --O GRCh38.exome.dict

# Index the reference for alignment with bwa
bwa index -a bwtsw GRCh38.exome.fa        
