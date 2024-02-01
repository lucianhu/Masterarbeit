# **DEVELOPING HUMAN REFERENCE PAN-EXOME FOR NEXT GENERATION SEQUENCING GENOMIC TESTING**
Master Thesis (Masterarbeit)

Submitted in Partial Fulfillment of the Degree Requirements Master of Science (M.Sc.) in Biotechnology

at the Department of Biotechnology, University of Applied Sciences Mannheim (Hochschule Mannheim), Germany.

by Quynh Nhu Nguyen

Under the Guidance of Dr. Phuc Loi Luu of “on-site Supervisor”, Grant & Innovation Center (GIC), University of Medicine & Pharmacy Ho Chi Minh City (UMP), Vietnam 

and Prof. Dr. Markus Gumbel, Hochschule Mannheim, Germany

Mannheim, February 2024 

This repository contains scripts used to reproduce our work.

# Workflow Description
The scripts are designed to be run in the order listed below:
- [Download raw data from the NCBI](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Download_raw_data)
- [Preparation GRCh38 reference genome](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Preparation_reference_genome)
- [Analysis nf-core/sarek pipeline to detect germline variants with the GRCh38 linear reference genome](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Linear_pipeline)
- [Preprocess VCFs and create consensus sequences](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Preprocess_VCFs_and_create%20_consensus%20_sequences)
- [Building pan-exome using nf-core/pangenome](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Building_pan-exome)
- [Up/Downstream pangenome WDL](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Pan-exome_pipeline)
- [Quantification and statistical analysis](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Statistics)

*Concise explanations of each command's purpose are embedded in comments throughout the script.*
