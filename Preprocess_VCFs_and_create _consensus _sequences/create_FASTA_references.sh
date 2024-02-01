#!/usr/bin/env bash
# From a human reference genome, create FASTA inputs for a Pangenome pipeline

# Set the path to the reference genome
REFERENCE='/media/lucianhu/29478145-8054-4fe3-900a-af5bbccd1109/Masterarbeit/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/hg38.non-exome.fa'

# Set the path to the interval BED file
INTERVAL_BED='/media/lucianhu/29478145-8054-4fe3-900a-af5bbccd1109/Masterarbeit/Homo_sapiens/GATK/GRCh38/Annotation/intervals/sorted_hg38-exome-interval-forcreateFASTA.bed'

# Enable debugging and handle errors
set -ex -o pipefail

# Create necessary folders/files
FASTA_DIR='/home/lucianhu/Pangenome/tumors/reference_sequences'
ERROR="${FASTA_DIR}/errors/references_Pan.log"
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
bedtools getfasta -fi "$REFERENCE" -bed "$INTERVAL_BED" -fo "REF.hg38.fa"

# Concatenate the sequences from all of FASTA entries into a single continuous sequence
awk '/^>/ {next} {printf("%s", $0)} END {printf("\n")}' "REF.hg38.fa" > "REF.hg38.PanSN.fa" && sed -i '1i>REF#chr' "REF.hg38.PanSN.fa"

# Compress a FASTA file
bgzip -@ 16 "REF.hg38.PanSN.fa" && samtools faidx "REF.hg38.PanSN.fa.gz"

echo "FASTA creation and indexing completed successfully."
