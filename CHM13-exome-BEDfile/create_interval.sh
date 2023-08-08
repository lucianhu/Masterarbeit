# Download GENCODE 38 genes of the T2T-CHM13v2.0 assembly.
aria2c -x8 https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3.gz -o CHM13-genes.gff3.gz

# Uncompress
gunzip *.gz

# Convert GTF/GXF file into bed file
agat_convert_sp_gff2bed.pl --gff CHM13-genes.gff3 --sub exon -o CHM13-exon.bed

# Sort BED file
awk -F'\t' 'BEGIN {OFS="\t"} {gsub(/^/, "chr", $1); print $1, $2, $3}' CHM13-exon.bed | grep -v -e '\t.\t' -e '^\.' -e '\t.$' > CHM13-interval.bed

sort -k1,1 -k2,2n < CHM13-interval.bed > sorted_CHM13-interval.bed