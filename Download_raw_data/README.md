## Step 1  
# SRA Sequences Download and Validation
Hereditary Breast and Ovarian Cancer Syndrome (HBOC) genomic data were acquired from twelve human NCBI Sequence Read Archive (SRA). These scripts are designed to facilitate the download and validation of SRA sequences. 
The process involves reading sample information from a provided CSV file, downloading corresponding SRA sequences, and then validating and renaming the obtained files.

Developmental and Epileptic Encephalopathies Syndrome (DEE) data, directly provided by the hospital, does not require running this script.
### Purpose:
- Download SRA sequences based on information in the provided CSV file.
- Rename and compress the downloaded files.
- Capture runtime metrics and reports any errors encountered during the process.
- Check if the paired-end FASTQ files have the expected number of reads.

### Prerequisites:
Ensure the following dependencies are installed before running the script:
- Install prefetch and fasterq-dump tools.
- Set the correct paths for the CSV_FILE and OUTPUT_DIR variables. 
CSV formatted for nf-core/sarek sample sheet usage.

### Usage:
```bash
bash Download_raw_data/1.downloadSRA.sh
bash Download_raw_data/2.validateSRA.sh
```






