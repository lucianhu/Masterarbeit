## Step 4

## Preprocess VCFs and create consensus sequences

### Purpose
The provided Bash scripts serve the purpose of merging variants from different tools into a single VCF file and creating FASTA inputs for a Pangenome pipeline. Additionally, there's a script to create FASTA sequences directly from a human reference genome for the Pangenome pipeline.

## Script 1: Merge Variants and Create FASTA Sequences

This script merges variants from different tools into a single VCF file and generates FASTA sequences for downstream Pangenome analysis. It processes input data specified in a CSV file containing a list of patient names.

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
This script creates FASTA inputs for a Pangenome pipeline directly from a human reference genome.

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
