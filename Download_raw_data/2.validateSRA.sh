#!/bin/env bash
# Validate SRA sequences

# Set the path to the destination folder
OUTPUT_DIR="path/to/raw_data/"
CSV_FILE="path/to/samplesheet_.csv"

# Read the list of files and sample names from the input CSV
tail -n +2 "$CSV_FILE" | while IFS=',' read -r PATIENT SAMPLE LANE FASTQ_1 FASTQ_2; do
    # Check if the renamed files exist
    if [ "$(find "$OUTPUT_DIR" -name "${SAMPLE}_*" | wc -l)" -eq 2 ]; then
        echo "Paired-end FASTQ ${SAMPLE} files have 2 reads each."
    else
        echo "Paired-end FASTQ ${SAMPLE} files do not have enough reads."
    fi
done

    


