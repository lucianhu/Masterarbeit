#!/usr/bin/env bash
# From a human reference genome, create FASTA inputs for a Pangenome pipeline

# Set the path to the reference genome
REFERENCE="path/to/GRCh38.exome.fa"

# Set the path to the interval BED file
INTERVAL_BED="/path/to/sorted_GRCh38-exome-interval.bed"

# Enable debugging and handle errors
set -ex -o pipefail

# Create necessary folders/files
FASTA_DIR="path/to/reference_sequences"
ERROR="${FASTA_DIR}/errors/references.log"
mkdir -p "${FASTA_DIR}/errors"

# Redirect both stdout and stderr to the log file
exec > "$ERROR" 2>&1

# Ensure we're working in the correct directory
cd "$FASTA_DIR" || exit 1

# Check if the reference genome and interval file exist
if [ ! -f "$REFERENCE" ]; then
    echo "Reference genome not found: $REFERENCE"
    exit 1
fi

if [ ! -f "$INTERVAL_BED" ]; then
    echo "Interval file not found: $INTERVAL_BED"
    exit 1
fi

# Extract target region
bedtools getfasta -fi "$REFERENCE" -bed "$INTERVAL_BED" -fo "REF.fa"

# Concatenate the sequences from all of FASTA entries into a single continuous sequence
awk '/^>/ {next} {printf("%s", $0)} END {printf("\n")}' "REF.fa" > "REF.PanSN.fa" && sed -i '1i>REF#chr' "REF.PanSN.fa"

# Compress a FASTA file
bgzip -@ 16 "REF.PanSN.fa" && samtools faidx "REF.PanSN.fa.gz"

echo "FASTA creation and indexing completed successfully."
