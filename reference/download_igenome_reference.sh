# Install AWS CLI
conda install -c conda-forge awscli=1.29.16

# Install aws-igenomes 
curl -fsSL https://ewels.github.io/AWS-iGenomes/aws-igenomes.sh > $HOME/DNA_softwares/aws-igenomes.sh

# Download aws-igenomes GATK_GRCh37
bash $HOME/DNA_softwares/aws-igenomes.sh -g Homo_sapiens -s GATK -b GRCh37 -o $HOME/Homo_sapiens/GATK/GRCh37 -q

# Download aws-igenomes GATK_GRCh38
bash $HOME/DNA_softwares/aws-igenomes.sh -g Homo_sapiens -s GATK -b GRCh38 -o $HOME/Homo_sapiens/GATK/GRCh38 -q 

# Download GATK_CHM13v2.0_Resource_Bundle
aws s3 --no-sign-request sync s3://human-pangenomics/T2T/CHM13/assemblies/variants/GATK_CHM13v2.0_Resource_Bundle/ $HOME/Homo_sapiens/GATK/CHM13/Annotation/GATKBundle

# Dowload CHM13 Fasta
aria2c -x8 https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -o $HOME/Homo_sapiens/GATK/CHM13/Sequence/WholeGenomeFasta/chm13v2.0.fa.gz

aria2c -x8 https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.gzi -o $HOME/Homo_sapiens/GATK/CHM13/Sequence/WholeGenomeFasta/chm13v2.0.fa.gz.gzi

