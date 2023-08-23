#!/bin/bash
# Complex variant comparison in the human genome reference HG19

if [ $# -ne 5 ]; then
    echo "Usage: $0 QUERY_DIR OUTPUT_DIR TRUTH BED STRATIFICATIONS"
    exit 1
fi

QUERY_DIR=$1

OUTPUT_DIR=$2

TRUTH=$3

BED=$4

STRATIFICATIONS=$5

REFERENCE='/mnt/sda/Masterarbeit/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/hg19.fa'

INTERVAL='/mnt/sda/Masterarbeit/Homo_sapiens/GATK/GRCh37/Annotation/intervals/sorted_hg19-exome-interval.bed'

VCFEVAL_PATH="$HOME/miniconda3/envs/hap.py/share/rtg-tools-3.12.1-0/rtg"

VCFEVAL_TEMPLATE='/mnt/sda/Masterarbeit/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/hg19_SDF'

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


