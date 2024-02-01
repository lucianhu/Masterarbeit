## Step 4

# Preprocess VCFs and create consensus sequences

## Script 1: Preprocess VCFs and Create FASTA Sequences

This script filters and phases variants, producing consensus sequences for constructing a pan-exome. It operates on input data provided in a CSV file containing a list of sample names. The generated FASTA files must adhere to the PanSN specification.

### Prerequisites

- Paths to relevant files and directories are correctly set.
- List of patients in CSV format `samplesheet.csv`.
- Reference genome `GRCh38.exome.fa`.
- Interval BED file `sorted_GRCh38-exome-interval.bed`.
- Segmentation duplications interval file `GRCh38_segdups_gt10kb.bed`.
- dbSNP file for variant validation `dbsnp_146.hg38.vcf.gz`.

### Usage

```bash
bash Preprocess_VCFs_and_create _consensus _sequences/2.execute_create_FASTA_samples.sh
```

## Script 2: Create FASTA from Reference Genome
This script creates a reference exome-FASTA input for constructing a pan-exome directly from a human reference genome. The FASTA file must conform to the PanSN specification. 

### Prerequisites
- Paths to the reference genome ("GRCh38.exome.fa") and interval BED file ("sorted_GRCh38-exome-interval.bed") are correctly set.

### Usage
```bash
bash Preprocess_VCFs_and_create _consensus _sequences/3.create_FASTA_references.sh
```

## Script 3: Concatenate FASTA Files
This script concatenates `*.PanSN.fa.gz` files into a single `pan-exome_FASTA.fa.gz` file, compresses it using bgzip, indexes it using `samtools faidx`, and calculates the length of each sequence.

### Prerequisites
- Input directory `samples_dir` containing `*.PanSN.fa.gz` files is available.

### Usage
```bash
bash Preprocess_VCFs_and_create _consensus _sequences/4.collect_PanSN_FASTAs.sh
```
