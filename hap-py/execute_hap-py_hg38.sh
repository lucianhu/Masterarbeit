#!/bin/bash
HAPPY_DIR="$HOME/masterarbeitlocal/Script/hap-py/hap-py_hg38.sh"

QUERY_DIR="$HOME/masterarbeitlocal/Benchmarks/HG38/HG001"

OUTPUT_DIR="${QUERY_DIR}/reports.Platinum-genomes"

TRUTH='/mnt/sda/Masterarbeit/Benchmarks/HG001/platinum-genomes/HG38/platinum-genomes_hg38_small-variants_NA12878.vcf.gz'

BED='/mnt/sda/Masterarbeit/Benchmarks/HG001/platinum-genomes/HG38/platinum-genomes_hg38_small-variants_NA12878.bed'

#OUTPUT_DIR="${QUERY_DIR}/report.NISTv4.2.1"

#TRUTH='/mnt/sda/Masterarbeit/Benchmarks/HG002/NISTv4.2.1/HG38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'

#BED='/mnt/sda/Masterarbeit/Benchmarks/HG002/NISTv4.2.1/HG38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'

#TRUTH='/mnt/sda/Masterarbeit/Benchmarks/HG002/CMRG_v1.00/HG38/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz'

#BED='/mnt/sda/Masterarbeit/Benchmarks/HG002/CMRG_v1.00/HG38/HG002_GRCh38_CMRG_smallvar_v1.00.bed'

STRATIFICATIONS='/mnt/sda/Masterarbeit/Benchmarks/genome-stratifications/HG001-GRCh38-all-stratifications.tsv'

bash ${HAPPY_DIR} ${QUERY_DIR} ${OUTPUT_DIR} ${TRUTH} ${BED} ${STRATIFICATIONS}


