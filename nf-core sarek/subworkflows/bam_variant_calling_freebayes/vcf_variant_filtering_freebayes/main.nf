include { GATK4_FILTERVARIANTFREEBAYES as FILTERVARIANTFREEBAYES } from '/home/lucianhu/.nextflow/assets/nf-core/sarek/modules/nf-core/freebayes/filtervariantfreebayes/main.nf'

workflow VCF_VARIANT_FILTERING_FREEBAYES {
    take:
    vcf             // channel: [ meta, vcf, tbi ]
    fasta           // channel: [ fasta ]
    fasta_fai       // channel: [ fasta_fai ]
    dict            // channel: [ dict ]
    
    main:
    versions = Channel.empty()

    freebayes_in = vcf.map{ meta, vcf, tbi -> [ meta, vcf, tbi ] }
    
    // Call the FILTERVARIANTFREEBAYES process
    FILTERVARIANTFREEBAYES(freebayes_in, fasta, fasta_fai, dict)
    
    // Get the filtered VCF output
    vcf = FILTERVARIANTFREEBAYES.out.vcf
        // Add variantcaller to the meta map and remove the num_intervals field
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'freebayes' ], vcf ] }
              
    emit:
    vcf
     
}
