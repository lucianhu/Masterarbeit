#!/bin/bash
# Convert VCF File

if [ $# -ne 2 ]; then
    echo "Usage: $0 QUERY_LIFTOVER LIFTED_OUT"
    exit 1
fi

## Activate the "rust+transanno" Conda environment
source activate rust

## Prepare a query FASTA, a reference FASTA, a chain file
QUERY_LIFTOVER=$1 # Should end with .gz
LIFTED_OUT=$2 # Should end with .gz

QUERY_FASTA='/mnt/sda/Masterarbeit/Homo_sapiens/CHM13/Sequence/WholeGenomeFasta/hs1.fa'
REFERENCE_FASTA='/mnt/sda/Masterarbeit/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/hg38.fa' 
CHAIN_FILE='/mnt/sda/Masterarbeit/Benchmarks/liftOver/hg38-chm13v2.over.chain.gz'

## Prepare VCF files
### Extract the directory path and extension from LIFTED_OUT
base_name=$(basename "${LIFTED_OUT}" .vcf.gz)
directory_path=$(dirname "${LIFTED_OUT}")
desired_part="${directory_path}/${base_name}"

REJECTED_OUT="${directory_path}/${base_name}.rejected.vcf.gz"
SORTED_OUT="${directory_path}/${base_name}.sorted.vcf.gz"
ERROR="${directory_path}/${base_name}.log"

## Run transanno to convert coordinates
cd $HOME/DNA_softwares/transanno

./target/release/transanno liftvcf -m --chain ${CHAIN_FILE} -o ${LIFTED_OUT} --query ${QUERY_FASTA} --reference ${REFERENCE_FASTA} --vcf ${QUERY_LIFTOVER} --fail ${REJECTED_OUT} > ${ERROR} 2>&1

## Sort output file
bcftools sort -O z -o ${SORTED_OUT} ${LIFTED_OUT} >> ${ERROR} 2>&1

bcftools index -t ${SORTED_OUT} >> ${ERROR} 2>&1

conda deactivate

