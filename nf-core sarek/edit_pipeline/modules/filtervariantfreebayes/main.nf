process GATK4_FILTERVARIANTFREEBAYES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path dict
    
    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.filtered.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: version

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK4_FILTERVARIANTFREEBAYES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    // Hard Filtering
    """
    # Extract SNPs and INDELS
    gatk --java-options "-Xmx${avail_mem}M" SelectVariants \\
        --reference $fasta \\
        --variant $vcf \\
        --select-type SNP \\
        --output ${prefix}.SNP.vcf.gz \\
        --tmp-dir .

    gatk --java-options "-Xmx${avail_mem}M" SelectVariants \\
        --reference $fasta \\
        --variant $vcf \\
        --select-type INDEL \\
        --output ${prefix}.INDEL.vcf.gz \\
        --tmp-dir .

    # Filter SNPs
    gatk --java-options "-Xmx${avail_mem}M" VariantFiltration \\
        --reference $fasta \\
        --variant ${prefix}.SNP.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD_filter" \\
        -filter "QUAL < 30.0" --filter-name "QUAL_filter" \\
        -filter "SOR > 3.0" --filter-name "SOR_filter" \\
        -filter "FS > 60.0" --filter-name "FS_filter" \\
        -filter "MQ < 40.0" --filter-name "MQ_filter" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum_filter" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_filter" \\
        -genotype-filter-expression "DP < 30" \\
        -genotype-filter-name "DP_filter" \\
        -genotype-filter-expression "GQ < 20" \\
        -genotype-filter-name "GQ_filter" \\
        --output ${prefix}.SNP.filtered.vcf.gz \\
        --tmp-dir .
     
    # Filter INDELs
    gatk --java-options "-Xmx${avail_mem}M" VariantFiltration \\
        --reference $fasta \\
        --variant ${prefix}.INDEL.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD_filter" \\
        -filter "QUAL < 30.0" --filter-name "QUAL_filter" \\
        -filter "SOR > 10.0" --filter-name "SOR_filter" \\
        -filter "FS > 200.0" --filter-name "FS_filter" \\
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum_filter" \\
        --output ${prefix}.INDEL.filtered.vcf.gz \\
        -genotype-filter-expression "DP < 30" \\
        -genotype-filter-name "DP_filter" \\
        -genotype-filter-expression "GQ < 20" \\
        -genotype-filter-name "GQ_filter" \\
        --tmp-dir .

    # Select Variants that "PASS" filters
    gatk --java-options "-Xmx${avail_mem}M" SelectVariants \\
        --exclude-filtered \\
        --variant ${prefix}.SNP.filtered.vcf.gz \\
        --output ${prefix}.SNP.ready.vcf.gz \\
        --tmp-dir .

    gatk --java-options "-Xmx${avail_mem}M" SelectVariants \\
        --exclude-filtered \\
        --variant ${prefix}.INDEL.filtered.vcf.gz \\
        --output ${prefix}.INDEL.ready.vcf.gz \\
        --tmp-dir .

    # Merge filtered SNPs and INDELs
    gatk --java-options "-Xmx${avail_mem}M" MergeVcfs \\
        --INPUT ${prefix}.SNP.ready.vcf.gz \\
        --INPUT ${prefix}.INDEL.ready.vcf.gz \\
        --OUTPUT ${prefix}.filtered.vcf.gz
                
    cat << END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """  
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

