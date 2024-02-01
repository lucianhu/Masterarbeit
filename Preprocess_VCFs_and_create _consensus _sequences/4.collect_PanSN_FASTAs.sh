#!/usr/bin/env bash

# Define the folder paths
INPUT_DIR="path/to/samples_dir"
DESTINATION_DIR="path/to/a_general_fasta_for_pangenome"

# Check if INPUT_DIR exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create the destination folder if it doesn't exist
mkdir -p "$DESTINATION_DIR"

# Navigate to the directory
cd "$INPUT_DIR"

# Use find to locate all '*.PanSN.fa' files under the 'samples' directory and its subdirectories,
# then concatenate them into 'BRCAs.fa'
find . -name '*.PanSN.fa.gz' -type f -exec zcat {} + > "$DESTINATION_DIR/pan-exome_FASTA.fa"

# Check if any '*.PanSN.fa' files were found
if [ ! -s "$DESTINATION_DIR/pan-exome_FASTA.fa" ]; then
    echo "No '*.PanSN.fa' files found to concatenate."
    exit 1
fi

# Compress *.fa using bgzip
bgzip -@ 16 "$DESTINATION_DIR/pan-exome_FASTA.fa"

# Check if bgzip compression was successful
if [ ! -s "$DESTINATION_DIR/pan-exome_FASTA.fa.gz" ]; then
    echo "Failed to compress 'pan-exome_FASTA.fa' using bgzip."
    exit 1
fi

# Index the compressed *.fa.gz file using samtools faidx
samtools faidx "$DESTINATION_DIR/pan-exome_FASTA.fa.gz"

# Check if indexing was successful
if [ ! -s "$DESTINATION_DIR/pan-exome_FASTA.fa.gz.fai" ]; then
    echo "Failed to index 'pan-exome_FASTA.fa.gz' using samtools."
    exit 1
fi

# Calculate the length of each sequence
zcat "$DESTINATION_DIR/pan-exome_FASTA.fa.gz"  | awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)} END {print seqlen}' > "$DESTINATION_DIR/length_of_each_sequence.txt"

echo "Concatenation, compression, and indexing completed successfully."
