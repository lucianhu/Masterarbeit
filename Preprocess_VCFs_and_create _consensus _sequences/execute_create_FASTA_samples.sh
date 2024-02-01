#!/usr/bin/env bash
# Merge variants from different tools into a single VCF file and create FASTA inputs for a Pangenome pipeline

# Set the path to the script that creates FASTA sequences from VCF files
CREATE_FASTA_BASH_DIR='/home/lucianhu/Mastersarbeit/Pangenome/FASTAs/create_FASTA_samples.sh'

# Set the path to the file containing the list of patient names
CSV_FILE="$HOME/Mastersarbeit/samplesheet_tumors.csv"

# Set the path to the reference genome
REFERENCE='/media/lucianhu/29478145-8054-4fe3-900a-af5bbccd1109/Masterarbeit/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/hg38.non-exome.fa'

# Set the path to the interval BED file
INTERVAL_BED='/media/lucianhu/29478145-8054-4fe3-900a-af5bbccd1109/Masterarbeit/Homo_sapiens/GATK/GRCh38/Annotation/intervals/sorted_hg38-exome-interval-forcreateFASTA.bed'

# Set the path to the region containing large segmental duplications 
SEGDUPS_INTERVAL='/media/lucianhu/29478145-8054-4fe3-900a-af5bbccd1109/Masterarbeit/Homo_sapiens/GATK/GRCh38/Annotation/intervals/GRCh38_segdups_gt10kb.forcreateFASTA.bed'

# Set the file path for the DBSNP file used for variant validation.
DBSNP_FILE='/media/lucianhu/29478145-8054-4fe3-900a-af5bbccd1109/Masterarbeit/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz'

# Check if the list of patients file exists
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: List of patients file not found: $CSV_FILE"
    exit 1
fi

# Loop through the list of patient names
tail -n +2 "$CSV_FILE" | while IFS=',' read -r PATIENT SAMPLE LANE FASTQ_1 FASTQ_2; do
    SAMPLE_NAME="$SAMPLE"
    VCF_FILE="/home/lucianhu/Masterarbeit/Linear_reference_genome/nf-core_sarek/execute/tumors/variant_calling/deepvariant/${SAMPLE_NAME}/${SAMPLE_NAME}.deepvariant.vcf.gz"   
    OUTPUT_FASTA_DIR="/home/lucianhu/Pangenome/tumors/${SAMPLE_NAME}/FASTAs"
    
    # Check if the VCF files exist
    if [ ! -f "$VCF_FILE" ]; then
        echo "Error: VCF file not found for $SAMPLE_NAME"
        continue  # Skip to the next patient
    fi

    # Check if the FASTA directory exists, create it if necessary
    mkdir -p "$OUTPUT_FASTA_DIR"

    echo "Processing $SAMPLE_NAME: Creating FASTA sequence..."

    # Execute the script to create FASTA sequences
    bash "$CREATE_FASTA_BASH_DIR" "$SAMPLE_NAME" "$VCF_FILE" "$OUTPUT_FASTA_DIR" "$REFERENCE" "$INTERVAL_BED" "$SEGDUPS_INTERVAL" "$DBSNP_FILE"

    echo "FASTA sequence creation for $SAMPLE_NAME completed."
done