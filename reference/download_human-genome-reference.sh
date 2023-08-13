#!/bin/bash

# HG19
## Download and uncompress file
aria2c -x8 https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.full_analysis_set.fa.gz -o hg19-raw.fa.gz

mv hg19-raw.fa.gz $HOME/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/hg19-raw.fa.gz

cd $HOME/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/

gzip -d hg19-raw.fa.gz

## Fix to be comparable with GATK resource bundle
awk '/^>chr/{sub("^>chr", ">");}1' hg19-raw.fa > hg19.fa

rm -f hg19-raw.fa

## Index fasta
samtools faidx hg19.fa

## Create Sequence Dictionary
gatk CreateSequenceDictionary \
--R hg19.fa \
--O hg19.dict

## Index the reference genome
bwa index -a bwtsw hg19.fa

mv hg19.amb hg19.ann hg19.bwt hg19.pac hg19.sa $HOME/Homo_sapiens/GATK/GRCh37/Sequence/BWAIndex

# HG38
## Download and uncompress file
aria2c -x8 https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.fullAnalysisSet.chroms.tar.gz

mv hg38.fullAnalysisSet.chroms.tar.gz $HOME/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/

cd $HOME/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/

tar xvzf hg38.fullAnalysisSet.chroms.tar.gz

cd hg38.fullAnalysisSet.chroms

cat *.fa > ../hg38.fa

##  Remove unnecessary files, except hg38.fa
cd ../

rm -rf hg38.fullAnalysisSet.chroms

rm -f hg38.fullAnalysisSet.chroms.tar.gz

## Index fasta
samtools faidx hg38.fa

## Create Sequence Dictionary
gatk CreateSequenceDictionary \
--R hg38.fa \
--O hg38.dict

## Index the reference genome
bwa index -a bwtsw hg38.fa

mv hg38.amb hg38.ann hg38.bwt hg38.pac hg38.sa $HOME/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex


# T2T-CHM13v2.0
## Download and uncompress file
aria2c -x8 https://hgdownload-test.gi.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz

mv hs1.fa.gz $HOME/Homo_sapiens/CHM13/Sequence/WholeGenomeFasta/hs1.fa.gz

cd $HOME/Homo_sapiens/CHM13/Sequence/WholeGenomeFasta

gzip -d hs1.fa.gz

## Index fasta
samtools faidx hs1.fa

## Create Sequence Dictionary
gatk CreateSequenceDictionary \
--R hs1.fa \
--O hs1.dict

## Index the reference genome
bwa index -a bwtsw hs1.fa

mv hs1.fa.amb hs1.fa.ann hs1.fa.bwt hs1.fa.pac hs1.fa.sa $HOME/Homo_sapiens/CHM13/Sequence/BWAIndex
