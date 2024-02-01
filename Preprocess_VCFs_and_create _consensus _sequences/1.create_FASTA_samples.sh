#!/usr/bin/env bash
# Preprocess, phase VCF files, and create haplotype FASTA inputs for a Pangenome pipeline

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 SAMPLE_NAME VCF_FILE OUTPUT_FASTA_DIR REFERENCE INTERVAL_BED SEGDUPS_INTERVAL(BED) DBSNP_FILE(VCF)"
    echo "Preprocess VCFs and create consensus FASTA sequence inputs for a Pangenome pipeline"
    exit 1
fi

# Assign command line arguments to variables
SAMPLE_NAME=$1
VCF_FILE=$2
OUTPUT_FASTA_DIR=$3
REFERENCE=$4
INTERVAL_BED=$5
SEGDUPS_INTERVAL=$6
DBSNP_FILE=$7

# Enable debugging and handle errors
set -exo pipefail

# Create necessary folders/files
mkdir -p "${OUTPUT_FASTA_DIR}/error" "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing" "${OUTPUT_FASTA_DIR}/temps"
ERROR="${OUTPUT_FASTA_DIR}/error/${SAMPLE_NAME}.log"

# Redirect both stdout and stderr to the log file
exec >>"$ERROR" 2>&1

# Ensure we're working in the correct directory
cd "${OUTPUT_FASTA_DIR}" || exit 1

# Preprocess VCFs: remove variants in large segmental duplications
rtg vcffilter -i "$VCF_FILE" -o "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.vcf.gz" --exclude-bed "$SEGDUPS_INTERVAL" --max-alleles=2 --min-alleles=2 --min-quality=30 --remove-overlapping --min-genotype-quality=20 --min-read-depth=35 --snps-only

# Preprocess VCFs: validate variants
gatk ValidateVariants -R "${REFERENCE}" -V "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.vcf.gz" --dbsnp "${DBSNP_FILE}" 

# Annotate VCF file with AC/AN fields filled
rtg vcfannotate -i "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.vcf.gz" -o "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.AC.AN.vcf.gz" --fill-an-ac

# Process each chromosome
parallel -j 10 --halt soon,fail=1 '
    i={}; 
    SAMPLE_NAME='"$SAMPLE_NAME"'; 
    OUTPUT_FASTA_DIR='"$OUTPUT_FASTA_DIR"'; 
    IN="${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chr${i}.vcf.gz";
    bcftools view "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.AC.AN.vcf.gz" --regions "chr${i}" -o $IN -Oz && tabix -p vcf $IN;
    MAP="path/to/phasing_reference_panel/genetic_maps.b38/chr${i}.b38.gmap.gz";
    CUT=0.001;
    THREADS=4;
    REFERENCE_PANEL="path/to/phasing_reference_panel/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz";
    OUT="${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chr${i}.phased.bcf";
    LOG="${OUTPUT_FASTA_DIR}/error/${SAMPLE_NAME}.chr${i}.phased.log";     
    SHAPEIT5_phase_common --input $IN --reference $REFERENCE_PANEL --map $MAP --region chr${i} --filter-maf $CUT --output $OUT --thread $THREADS --log $LOG --progress && bcftools index -f $OUT --threads $THREADS;
    bcftools convert -Oz -o "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chr${i}.phased.vcf.gz" $OUT;
    tabix -p vcf "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chr${i}.phased.vcf.gz";
    rm -f $IN "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chr${i}.vcf.gz.tbi" $OUT "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chr${i}.phased.bcf.csi";
' ::: {1..22} X

# Process chromosome Y separately
bcftools view "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.AC.AN.vcf.gz" --regions "chrY" -Oz -o "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chrY.vcf.gz"
tabix -p vcf "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.chrY.vcf.gz"

# Process homozygous variants separately
bcftools view "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.AC.AN.vcf.gz" -i 'GT="1/1"' -Oz -o "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.hom.vcf.gz"
tabix -p vcf "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/${SAMPLE_NAME}.no_dups.AC.AN.hom.vcf.gz"

# Merge phased VCFs
rtg vcfmerge --force-merge-all -o "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.merge.phased.vcf.gz" "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing/"*.vcf.gz

# Apply variants to create consensus sequence
bcftools consensus -H 1 -f "${REFERENCE}" "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.merge.phased.vcf.gz" > "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype1.full_hg38.fa" && samtools faidx "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype1.full_hg38.fa"

bcftools consensus -H 2 -f "${REFERENCE}" "${OUTPUT_FASTA_DIR}/preprocess_VCFs/${SAMPLE_NAME}.no_dups.merge.phased.vcf.gz" > "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype2.full_hg38.fa" && samtools faidx "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype2.full_hg38.fa"

# Extract target region
bedtools getfasta -fi "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype1.full_hg38.fa" -bed "$INTERVAL_BED" -fo "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype1.target_region.fa"

bedtools getfasta -fi "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype2.full_hg38.fa" -bed "$INTERVAL_BED" -fo "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype2.target_region.fa"

# Create a continuous sequence for haplotype1
awk '/^>/ {next} {printf("%s", $0)} END {printf("\n")}' "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype1.target_region.fa" > "${SAMPLE_NAME}.haplotype1.target_region.PanSN.fa" && sed -i "1i>${SAMPLE_NAME}_1#chr" "${SAMPLE_NAME}.haplotype1.target_region.PanSN.fa"

# Create a continuous sequence for haplotype2
awk '/^>/ {next} {printf("%s", $0)} END {printf("\n")}' "${OUTPUT_FASTA_DIR}/temps/${SAMPLE_NAME}.haplotype2.target_region.fa" > "${SAMPLE_NAME}.haplotype2.target_region.PanSN.fa" && sed -i "1i>${SAMPLE_NAME}_2#chr" "${SAMPLE_NAME}.haplotype2.target_region.PanSN.fa"

# Compress the created FASTA files
bgzip -@ 16 "${SAMPLE_NAME}.haplotype1.target_region.PanSN.fa" && samtools faidx "${SAMPLE_NAME}.haplotype1.target_region.PanSN.fa.gz"

bgzip -@ 16 "${SAMPLE_NAME}.haplotype2.target_region.PanSN.fa" && samtools faidx "${SAMPLE_NAME}.haplotype2.target_region.PanSN.fa.gz"

# Optionally, you can add messages to indicate the completion of each step
echo "FASTA sequence creation for ${SAMPLE_NAME} completed."

# Remove unnecessary files
rm -rf "${OUTPUT_FASTA_DIR}/preprocess_VCFs/phasing" "${OUTPUT_FASTA_DIR}/temps"
