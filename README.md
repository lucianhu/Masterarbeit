# **OPTIMIZATION OF WHOLE EXOME SEQUENCE PIPELINE: HUMAN REFERENCE PAN/GENOME**
Master Thesis (Masterarbeit)

Submitted in Partial Fulfillment of the Degree Requirements Master of Science (M.Sc.) in Biotechnology

at the Department of Biotechnology, University of Applied Sciences Mannheim (Hochschule Mannheim), Germany.

by Quynh Nhu Nguyen

Under the Guidance of Dr. Phuc Loi Luu of “on-site Supervisor” and Prof. Dr. Markus Gumbel, Hochschule Mannheim.

Mannheim, February 2024 

This repository contains scripts used to reproduce our work.

# Workflow Description
The scripts are designed to be run in the order listed below:
- Download raw data from the NCBI.
- Preparation GRCh38 reference genome
- Analysis nf-core/sarek pipeline to detect germline variants with the GRCh38 linear reference genome
- Preprocess VCFs and create consensus sequences
- Building WES pangenome graphs using nf-core/pangenome
- Up/Downstream pangenome WDL

*Concise explanations of each command's purpose are embedded in comments throughout the script.*
