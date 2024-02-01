## Step 3
# Analysis nf-core/sarek pipeline to detect germline variants with the GRCh38 linear reference genome

The purpose of this command is to run the `nf-core/sarek` workflow using specific configurations. It leverages version 3.4.0 of the workflow, utilizes the Docker profile for containerization, sets the working directory to `path/to/work`, and provides additional parameters from the specified YAML file `path/to/nf-params_Linear.yaml`.

### Prerequisites

Ensure the following prerequisites are met before executing the command:

- Nextflow is installed on your system.
- Docker is installed and available for containerization.
- The nf-core/sarek workflow (version 3.4.0) is accessible.

### Usage
```bash
nextflow run nf-core/sarek -r 3.4.0 -profile docker -work-dir "path/to/work" -params-file "path/to/nf-params_Linear.yaml"
```
with params.yaml containing
```
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh38'
```
and samplesheet.csv
```
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```
Our `params.yaml` configuration is accessible in [Linear_pipeline/nf-params_Linear.yaml](https://github.com/LuciaNhuNguyen/Masterarbeit/blob/main/Linear_pipeline/nf-params_Linear.yaml), while our samplesheets are located at [main/samplesheet_DEE.csv](https://github.com/LuciaNhuNguyen/Masterarbeit/blob/main/samplesheet_DEE.csv) and [main/samplesheet_tumors.csv](https://github.com/LuciaNhuNguyen/Masterarbeit/blob/main/samplesheet_tumors.csv). 

For any additional information, particularly concerning parameters related to the pipeline, please feel free to contact the [nf-core/sarek Team](https://nf-co.re/sarek/3.4.0/parameters#snpeff_cache).

### Note
The nf-core/sarek pipeline is optimized for AWS cloud infrastructure with iGenomes references. Local execution may encounter performance challenges; therefore, it is essential to download and locally store iGenomes references. Ensure that the link in your param.yaml is updated accordingly. This not only enhances performance but also allows for flexibility in selecting a different target reference genome, as illustrated in our case, where we can switch to `GRCh38.exome.fa`
```bash
# Install aws-igenomes 
curl -fsSL https://ewels.github.io/AWS-iGenomes/aws-igenomes.sh > aws-igenomes.sh

# Download aws-igenomes GATK_GRCh38
bash aws-igenomes.sh -g Homo_sapiens -s GATK -b GRCh38 -o "path/to/Homo_sapiens/GATK/GRCh38" -q
```
