#!/bin/bash
HAPPY_DIR='/home/lucianhu/masterarbeitlocal/Script/hap-py/hap-py_chm13.sh'

QUERY_DIR='/home/lucianhu/masterarbeitlocal/Benchmarks/CHM13/HG001'

OUTPUT_DIR="${QUERY_DIR}/reports.Platinum-genomes"

TRUTH='/mnt/sda/Masterarbeit/Benchmarks/HG001/platinum-genomes/CHM13/platinum-genomes_hg38-chm13_small-variants_NA12878.sorted.vcf.gz'

BED='/mnt/sda/Masterarbeit/Benchmarks/HG001/platinum-genomes/CHM13/platinum-genomes_hg38-chm13_small-variants_NA12878.bed'

#TRUTH='/mnt/sda/Masterarbeit/Benchmarks/HG002/CMRG_v1.00/CHM13/HG002_hg38-chm13_CMRG_smallvar_v1.00.sorted.vcf.gz'

#BED='/mnt/sda/Masterarbeit/Benchmarks/HG002/CMRG_v1.00/CHM13/HG002_hg38-chm13_CMRG_smallvar_v1.00.bed'

bash ${HAPPY_DIR} ${QUERY_DIR} ${OUTPUT_DIR} ${TRUTH} ${BED}


#TRUTH='/mnt/sda/Masterarbeit/Benchmarks/HG002/NISTv4.2.1/CHM13/HG002_hg38-chm13_NISTv4.2.1.transanno.sorted.vcf.gz'
#BED='/mnt/sda/Masterarbeit/Benchmarks/HG002/NISTv4.2.1/CHM13/HG002_hg38-chm13_NISTv4.2.1.benchmark.bed'