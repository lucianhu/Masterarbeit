process SNPEFF_SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf)
    val db
    path cache

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    path "*.csv", emit: report
    path "*.html", emit: summary_html
    path "*.genes.txt", emit: genes_txt
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    def snpEffJarPath = '$HOME/DNA_softwares/snpEff/snpEff.jar' // Replace with the actual path to snpEff.jar on your local machine
    def snpSiftJarPath = '$HOME/DNA_softwares/snpEff/SnpSift.jar' // Replace with the actual path to snpSift.jar on your local machine
    def clinvar = "$HOME/DNA_softwares/snpEff/data/${db}/${db}.clinvar.vcf.gz" // Correct the path to ClinVar.vcf.gz

    """
    java -Xmx${avail_mem}M -jar $snpEffJarPath $db $vcf -csvStats ${prefix}.csv > ${prefix}.ann.vcf

    gatk VariantsToTable -V ${prefix}.ann.vcf -F CHROM -F POS -F TYPE -F ID -F ANN -F LOF -F NMD -GF AD -GF DP -GF GQ -GF GT -O ${prefix}.ann.vcf.csv
    
    java -Xmx${avail_mem}M -jar $snpSiftJarPath annotate $clinvar $vcf > ${prefix}.ann.clinvar.vcf

    cp ${prefix}.ann.clinvar.vcf ${prefix}.ann.clinvar.vcf.csv

    gatk VariantsToTable -V ${prefix}.ann.clinvar.vcf -F CHROM -F POS -F TYPE -F ID -F ALLELEID -F CLNDN -F CLNSIG -F CLNSIGCONF -F CLNSIGINCL -F CLNVC -F GENEINFO -GF AD -GF GQ -GF GT -O ${prefix}.ann.clinvar.csv

    snpEff_version="NFCORE_SAREK:SAREK:VCF_ANNOTATE_ALL:VCF_ANNOTATE_SNPEFF:SNPEFF_SNPEFF \$(java -jar \"$snpEffJarPath\" -version | cut -f 2 -d ' ')" > versions.yml
   
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf
    
    snpEff_version="NFCORE_SAREK:SAREK:VCF_ANNOTATE_ALL:VCF_ANNOTATE_SNPEFF:SNPEFF_SNPEFF \$(java -jar \"$snpEffJarPath\" -version | cut -f 2 -d ' ')" > versions.yml 
    """
}

