#!/bin/bash
# Complex variant comparison in the human reference T2T-CHM13v2.0

if [ $# -ne 4 ]; then
    echo "Usage: $0 QUERY_DIR OUTPUT_DIR TRUTH BED"
    exit 1
fi

QUERY_DIR=$1

OUTPUT_DIR=$2

TRUTH=$3

BED=$4

REFERENCE='/mnt/sda/Masterarbeit/Homo_sapiens/CHM13/Sequence/WholeGenomeFasta/hs1.fa'

INTERVAL='/mnt/sda/Masterarbeit/Homo_sapiens/CHM13/Annotation/intervals/sorted_CHM13-exome-interval.bed'

STRATIFICATIONS='/mnt/sda/Masterarbeit/Benchmarks/genome-stratifications/CHM13_stratifications.tsv'

VCFEVAL_PATH='/home/lucianhu/miniconda3/envs/hap.py/share/rtg-tools-3.12.1-0/rtg'

VCFEVAL_TEMPLATE='/mnt/sda/Masterarbeit/Homo_sapiens/CHM13/Sequence/WholeGenomeFasta/hs1_SDF'

## Create the destination folder if it doesn't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p "${OUTPUT_DIR}/error"

## Activate the hap.py environment
source activate hap.py

## Run hap.py

for VCF_FILE in ${QUERY_DIR}/*.vcf.gz; do
    QUERY=$(readlink -f ${VCF_FILE})
    BASE_NAME=$(basename ${VCF_FILE})
    BASE_NAME_NOEXT="${BASE_NAME%.vcf.gz}"
    
    cd ${OUTPUT_DIR}
        
    ERROR="${OUTPUT_DIR}/error/${BASE_NAME_NOEXT}.log"
    
    hap.py ${TRUTH} ${QUERY} -r ${REFERENCE} -f ${BED} --threads 4 -T ${INTERVAL} -o "${BASE_NAME_NOEXT}.reports" \
    --preprocess-truth \
    --usefiltered-truth \
    --stratification ${STRATIFICATIONS} \
    --roc QUAL --roc-filter LowQual \
    --no-decompose --no-leftshift --pass-only \
    --engine vcfeval --engine-vcfeval-path ${VCFEVAL_PATH} --engine-vcfeval-template ${VCFEVAL_TEMPLATE} > ${ERROR} 2>&1
done

## Deactivate the environment

conda deactivate


