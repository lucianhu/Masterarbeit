# Go to SnpEff's install dir
cd $HOME/DNA_softwares/snpEff

# Create database dir
mkdir data/CHM13.2
cd data/CHM13.2

# Download genome
aria2c -x8 https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -o sequences.fa.gz
gzip -d sequences.fa.gz
samtools faidx sequences.fa

# Download annotated genes
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz -o genes.gtf.gz

# Download proteins
# This is used for:
#   - "Rare Amino Acid" annotations
#   - Sanity check (checking protein predicted from DNA sequences match 'real' proteins)
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-pep.fa.gz -o protein.fa.gz

# Download CDSs
# Note: This is used as "sanity check" (checking that CDSs predicted from gene sequences match 'real' CDSs)
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-cds.fa.gz -o cds.fa.gz

# Uncompress
gunzip *.gz

# Download ClinVar
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-clinvar.vcf.gz -o chm13v2.0_ClinVar.vcf.gz
gzip -d chm13v2.0_ClinVar.vcf.gz
bgzip chm13v2.0_ClinVar.vcf && tabix -p vcf chm13v2.0_ClinVar.vcf.gz

# Uncompress:
# Why do we need to uncompress?
# Because ENSEMBL compresses files using a block compress gzip which is not compatible with Java's library Gunzip

# Edit snpEff.config file
#
# WARNING! You must do this yourself. Just copying and pasting this into a terminal won't work.
#
# Add lines:
# T2T-CHM13v2.0 release from ENSEMBL
# CHM13.2.genome : Human genome (Homo sapiens) T2T-CHM13v2.0 using transcripts
# CHM13.2.reference: https://projects.ensembl.org/hprc/

# Now we are ready to build the database
java -Xmx20g -jar $HOME/DNA_softwares/snpEff/snpEff.jar build -v CHM13.2 2>&1 | tee CHM13.2.build
