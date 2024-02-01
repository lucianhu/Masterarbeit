#!/bin/bash

destination="/mnt/sda/Masterarbeit/raw_data/"

# Read the list of files and sample names from the input tsv
while IFS=$'\t' read -r sra name_of_sample; do
    
    # Find the files with the pattern SRRxxxxxxx_1.fastq.gz and rename them
    find "$destination" -name "${sra}_1.fastq.gz" -exec mv {} "$destination/${name_of_sample}_1.fastq.gz" \;

    # Find the files with the pattern SRRxxxxxxx_2.fastq.gz and rename them
    find "$destination" -name "${sra}_2.fastq.gz" -exec mv {} "$destination/${name_of_sample}_2.fastq.gz" \;

    # Check if the renamed files exist
    if [ "$(find "$destination" -name "${name_of_sample}_*" | wc -l)" -eq 2 ]; then
        echo "Paired-end FASTQ ${name_of_sample} files have enough 2 reads."
    else
        echo "Paired-end FASTQ ${name_of_sample} files have not enough."
    fi
    
done < list_SRA_raw_data.tsv


