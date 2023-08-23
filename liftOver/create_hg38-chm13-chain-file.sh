#!/bin/bash
# Create chain file from minimap2 result

## Activate the "rust+transanno" Conda environment
source activate rust

## Prepare files
QUERY_FASTA='/mnt/sda/Masterarbeit/Homo_sapiens/CHM13/Sequence/WholeGenomeFasta/hs1.fa'

REFERENCE_FASTA='/mnt/sda/Masterarbeit/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/hg38.fa' 

PAF_FILE='/mnt/sda/Masterarbeit/Benchmarks/liftOver/hg38-chm13.transanno.paf'

CHAIN_FILE='/mnt/sda/Masterarbeit/Benchmarks/liftOver/hg38-chm13.transanno.chain'

ERROR='/mnt/sda/Masterarbeit/Benchmarks/liftOver/error_chainfile.transanno.log'

## Run minimap2
minimap2 -cx asm5 --cs ${QUERY_FASTA} ${REFERENCE_FASTA} > ${PAF_FILE} 2> ${ERROR}

## Run transanno to create chain file and log errors
cd $HOME/DNA_softwares/transanno

./target/release/transanno minimap2chain ${PAF_FILE} --output ${CHAIN_FILE} >> ${ERROR} 2>&1

conda deactivate