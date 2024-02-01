# Step 5

## Building pan-exome using nf-core/pangenome

### Purpose
The purpose of the provided commands is as follows:
- **Nextflow Execution:**
  Initiates the `nf-core/pangenome` workflow, version 1.0.0, using Docker as the profile, with a specified working directory and parameter file.
- **Autoindexing with `vg autoindex`:**
  Generates an autoindexed representation of the pangenome using the Giraffe workflow, producing relevant output files.
- **Extracting Paths with `vg paths`:**
  Extracts GRCh38 reference path from the generated pangenome, specifically the `REF` path, and saves the list in the file `Pan-exome.reference_path.txt`.

# Prerequisites

Before executing these commands, ensure the following prerequisites are met:
- `ǹextflow` is installed on your system.
- `docker` is available for containerization.
- The necessary data files for the pangenome workflow are accessible.
- The `vg` tool (VariationGraph toolkit) is installed.

# Usage
Execute the commands in sequence:

Run the `nf-core/pangenome` workflow. Input data: `nf-params-Pan.yaml` containing the absolute path of a bgzip-compressed FASTA file containing all assemblies, and the FASTA file must adhere to the PanSN-spec. 
Output data: `Pan-exome.fa.gz.gfaffix.unchop.Ygs.view.gfa`.
```bash
nextflow run nf-core/pangenome -r 1.0.0 -profile docker -work-dir "path/to/work" -params-file "path/to/nf-params-Pan.yaml"
```

Autoindex with `vg autoindex`. Input data: `Pan-exome.fa.gz.gfaffix.unchop.Ygs.view.gfa`. Output data: `Pan-exome.giraffe.gbz`, `Pan-exome.dist`, `Pan-exome.min`.
```bash
vg autoindex --prefix Pan-exome --workflow giraffe --gfa Pan-exome.fa.gz.gfaffix.unchop.Ygs.view.gfa
```
Extract paths with `vg paths`. Input data: `Pan-exome.giraffe.gbz`. Output data: `Pan-exome.reference_path.txt`.
```bash
vg paths -S "REF" --list -x Pan-exome.giraffe.gbz > Pan-exome.reference_path.txt
````
