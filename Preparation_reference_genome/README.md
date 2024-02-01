## Step 2  
# Preparation GRCh38 reference genome

The purpose of this script is to automate the processing of the human reference genome (GRCh38) for genomic analysis and alignment pipelines. It performs essential tasks such as downloading the genome files, preparing the reference genome, indexing, extracting introns, creating a masked genome, and generating a sequence dictionary.

# Prerequisites

Ensure the following dependencies are installed before running the script:

- `wget` for downloading files
- `gunzip` for unzipping files
- `samtools` for indexing and handling sequence files
- `bedtools` for genomic interval operations
- `bwa` for creating an index for alignment
- `gatk` (Genome Analysis Toolkit) for creating a sequence dictionary

# Usage

Execute the script using the following command:

```bash
bash Preparation_reference_genome/1.make_reference.sh
```
