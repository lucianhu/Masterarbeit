#!/bin/bash

# HG19
## Download Illumina Exome 2.0 Plus Panel HG19 BED File
aria2c -x8 https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg19_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED

## Sort BED file
sort -k1,1 -k2,2n < hg19_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED > sorted_hg19-raw-exome-interval.bed

## Fix to be comparable with GATK resource bundle
awk '{ sub(/^chr/, "", $1); print $1 "\t" $2 "\t" $3 "\t" $4 }' sorted_hg19-raw-exome-interval.bed > sorted_hg19-exome-interval.bed

## Remove unnecessary files
rm -f hg19_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED
rm -f sorted_hg19-raw-exome-interval.bed

# HG38
## Download Illumina Exome 2.0 Plus Panel HG19 BED File
aria2c -x8 https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED

## Sort BED file
sort -k1,1 -k2,2n < hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED > sorted_hg38-exome-interval.bed

## Remove unnecessary files
rm -f hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED

# T2T-CHM13v2.0
## Download GENCODE 38 genes of the T2T-CHM13v2.0 assembly.
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3.gz -o CHM13-genes.gff3.gz

## Uncompress
gunzip *.gz

## Convert GTF/GXF file into bed file
agat_convert_sp_gff2bed.pl --gff CHM13-genes.gff3 --sub exon -o CHM13-exome.bed

## Fix to be comparable with UCSC resource
awk -F'\t' 'BEGIN {OFS="\t"} {gsub(/^/, "chr", $1); print $1, $2, $3}' CHM13-exome.bed | grep -v -e '\t.\t' -e '^\.' -e '\t.$' > CHM13-interval-exome.bed

## Sort BED file
sort -k1,1 -k2,2n < CHM13-interval-exome.bed > sorted_CHM13-exome-interval.bed

## Remove unnecessary files
rm -f CHM13-genes.gff3
rm -f CHM13-exome.bed
rm -f CHM13-interval-exome.bed

