#!/bin/bash
LIFTOVERVCF='/home/lucianhu/masterarbeitlocal/Script/liftOver/liftOver_vcf_hg38-chm13.sh'

#QUERY_LIFTOVER='/mnt/sda/Masterarbeit/Benchmarks/HG002/CMRG_v1.00/HG38/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz'
#LIFTED_OUT='/mnt/sda/Masterarbeit/Benchmarks/HG002/CMRG_v1.00/CHM13/HG002_hg38-chm13_CMRG_smallvar_v1.00.vcf.gz'

QUERY_LIFTOVER='/mnt/sda/Masterarbeit/Benchmarks/HG001/platinum-genomes/HG38/platinum-genomes_hg38_small-variants_NA12878.vcf.gz'
LIFTED_OUT='/mnt/sda/Masterarbeit/Benchmarks/HG001/platinum-genomes/CHM13/platinum-genomes_hg38-chm13_small-variants_NA12878.vcf.gz'

bash ${LIFTOVERVCF} ${QUERY_LIFTOVER} ${LIFTED_OUT}

