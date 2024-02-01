version 1.0

################################################################
# PREPARE DATA #
################################################################

task EXTRACTSUBSETPATHNAMES {
    input {
        File in_gbz_file
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G")) + 20
        Int in_extract_mem = 120
    }

    command {
        set -eux -o pipefail

        vg gbwt -CL -Z ${in_gbz_file} | sort > path_list.txt

        grep -v _decoy path_list.txt | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM > path_list.sub.txt
    }
    output {
        File output_path_list_file = "path_list.sub.txt"
    }
    runtime {
        preemptible: 2
        time: 30
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.53.0"
    }
}

task EXTRACTREFERENCE {
    input {
        File in_gbz_file
        File in_path_list_file
        String in_prefix_to_strip = ""
        Int in_extract_mem = 120
        Int in_extract_disk = 2 * round(size(in_gbz_file, "G")) + 20
    }

    command {
        set -eux -o pipefail

        # Subset to just the paths we care about (may be the whole file) so we
        # get a good dict with just those paths later
        vg paths \
           --extract-fasta \
           -p ${in_path_list_file} \
           --xg ${in_gbz_file} > ref.fa
        
        if [ ~{in_prefix_to_strip} != "" ]
        then
            mv ref.fa ref.prefix.fa
            sed -e "s/>~{in_prefix_to_strip}/>/g" ref.prefix.fa > ref.fa
        fi
    }
    output {
        File reference_file = "ref.fa"
    }
    runtime {
        preemptible: 2
        time: 100
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.53.0"
    }
}

task INDEXREFERENCE {
    input {
        File in_reference_file
        Int in_mem = 4
        Int in_disk = round(2 * size(in_reference_file, "G")) + 10
        Int timeMinutes = 1 + ceil(size(in_reference_file, "G"))
    }

    command <<<
        set -eux -o pipefail
      
        ln -s ~{in_reference_file} ref.fa
                
        # Index the subset reference
        samtools faidx ref.fa 
        
        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R=ref.fa \
          O=ref.dict

    >>>
    output {
        File reference_index_file = "ref.fa.fai"
        File reference_dict_file = "ref.dict"
    }
    runtime {
        preemptible: 2
        time: timeMinutes
        cpu: 1
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task CONVERT_CRAM_TO_FASTQ {
    input {
	    File? in_cram_file
        File? in_ref_file
        File? in_ref_index_file
        Boolean in_paired_reads = true
	    Int in_cores
        Int timeMinutes = 1 + ceil(size(in_cram_file, "G"))
        Int in_disk = round(5 * size(in_cram_file, "G")) + 50
    }
    Int half_cores = in_cores / 2

    command <<<
    set -eux -o pipefail

    if [ ~{in_paired_reads} == true ]
    then
        samtools collate -@ ~{half_cores} --reference ~{in_ref_file} -Ouf ~{in_cram_file} | samtools fastq -@ ~{half_cores} -1 reads.R1.fastq.gz -2 reads.R2.fastq.gz -0 reads.o.fq.gz -s reads.s.fq.gz -c 1 -N -
    else
        samtools fastq -@ ~{in_cores} -o reads.R1.fastq.gz -c 1 --reference ~{in_ref_file} ~{in_cram_file}
    fi

    >>>
    output {
        File output_fastq_1_file = "reads.R1.fastq.gz"
        File? output_fastq_2_file = "reads.R2.fastq.gz"
    }
    runtime {
        preemptible: 2
        time: timeMinutes
        cpu: in_cores
        memory: "50 GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    }    
}

################################################################
# SEQUENCING QUALITY CONTROL #
################################################################

task FASTQC_SINGLE_END {
    input {
        String in_sample_name
        File in_read_1_fastq
        Int in_mem = 4
        Int thread_count = 2
        Int in_disk = round(3 * size(in_read_1_fastq, "G")) + 50 
        Int timeMinutes = 1 + ceil(size(in_read_1_fastq, "G")) * 4   
    }

    command <<<
        set -eux -o pipefail
        
        [ ! -f ~{in_sample_name}_1.fastq.gz ] && ln -sf ~{in_read_1_fastq} ~{in_sample_name}_1.gz
                
        fastqc --quiet --threads ~{thread_count} ~{in_sample_name}_1.gz

    >>>
    output {
        File output_fastqc_html_1 = "~{in_sample_name}_1_fastqc.html"
        File output_fastqc_report_zip_1 = "~{in_sample_name}_1_fastqc.zip"            
    }
    runtime {
        time: timeMinutes
        cpu: thread_count
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    }
}

task FASTQC_PAIRED_READS {
    input {
        String in_sample_name
        File in_read_1_fastq
        File in_read_2_fastq
        Int in_mem = 4
        Int thread_count = 2
        Int in_disk = round(3 * size(in_read_1_fastq, "G")) + 50 
        Int timeMinutes = 1 + ceil(size(in_read_1_fastq, "G")) * 4   
    }

    command <<<
        set -eux -o pipefail

        [ ! -f ~{in_sample_name}_1.fastq.gz ] && ln -sf ~{in_read_1_fastq} ~{in_sample_name}_1.gz
        [ ! -f ~{in_sample_name}_2.fastq.gz ] && ln -sf ~{in_read_2_fastq} ~{in_sample_name}_2.gz
        
        fastqc --quiet --threads ~{thread_count} ~{in_sample_name}_1.gz ~{in_sample_name}_2.gz

    >>>
    output {
        File output_fastqc_html_1 = "~{in_sample_name}_1_fastqc.html"
        File output_fastqc_report_zip_1 = "~{in_sample_name}_1_fastqc.zip"
        File output_fastqc_html_2 = "~{in_sample_name}_2_fastqc.html"
        File output_fastqc_report_zip_2 = "~{in_sample_name}_2_fastqc.zip"       
    }
    runtime {
        time: timeMinutes
        cpu: thread_count
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    }
}

task FASTP_SINGLE_END {
    input {
        String in_sample_name
        File in_read_1_fastq
        Int in_call_mem
        Int thread_count = 12
        Int in_disk = round(3 * size(in_read_1_fastq, "G")) + 50
        Int timeMinutes = 1 + ceil(size(in_read_1_fastq, "G")  * 6.0 / thread_count)       
    }

    command <<<
        set -eux -o pipefail

        [ ! -f ~{in_sample_name}_1.fastq.gz ] && ln -sf ~{in_read_1_fastq} read_1.fastq.gz
                
        fastp \
          --in1 read_1.fastq.gz \
          --out1 ~{in_sample_name}_1.fastp.fastq.gz \
          --json ~{in_sample_name}.fastp.json \
          --html ~{in_sample_name}.fastp.html \
          --overrepresentation_analysis --correction --detect_adapter_for_pe \
          2> ~{in_sample_name}.fastp.log
    
    >>>
    output {
        File output_fastp_1_file = "~{in_sample_name}_1.fastp.fastq.gz"
        File log = "~{in_sample_name}.fastp.log"
        File json = "~{in_sample_name}.fastp.json"
        File html = "~{in_sample_name}.fastp.html"      
    }
    runtime {
        time: timeMinutes
        cpu: thread_count
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
    }
}

task FASTP_PAIRED_READS {
    input {
        String in_sample_name
        File in_read_1_fastq
        File in_read_2_fastq
        Int in_call_mem
        Int thread_count = 12
        Int in_disk = round(3* size(in_read_1_fastq, "G")) + 50 
        Int timeMinutes = 1 + ceil(size(in_read_1_fastq, "G")  * 6.0 / thread_count)    
    }

    command <<<
        set -eux -o pipefail

        [ ! -f ~{in_sample_name}_1.fastq.gz ] && ln -sf ~{in_read_1_fastq} read_1.fastq.gz
        [ ! -f ~{in_sample_name}_2.fastq.gz ] && ln -sf ~{in_read_2_fastq} read_2.fastq.gz
        
        fastp \
          --in1 read_1.fastq.gz \
          --in2 read_2.fastq.gz \
          --out1 ~{in_sample_name}_1.fastp.fastq.gz \
          --out2 ~{in_sample_name}_2.fastp.fastq.gz \
          --json ~{in_sample_name}.fastp.json \
          --html ~{in_sample_name}.fastp.html \
          --overrepresentation_analysis --correction --detect_adapter_for_pe \
          2> ~{in_sample_name}.fastp.log
    
    >>>
    output {
        File output_fastp_1_file = "~{in_sample_name}_1.fastp.fastq.gz"
        File output_fastp_2_file = "~{in_sample_name}_2.fastp.fastq.gz"
        File log = "~{in_sample_name}.fastp.log"
        File json = "~{in_sample_name}.fastp.json"
        File html = "~{in_sample_name}.fastp.html"      
    }
    runtime {
        time: timeMinutes
        cpu: thread_count
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
    }
}

task SPLIT_READS {
    input {
        File in_read_file
        String in_pair_id
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk = round(5 * size(in_read_file, "G")) + 20
    }

    command <<<
        set -eux -o pipefail

        CHUNK_LINES=$(( ~{in_reads_per_chunk} * 4 ))
        gzip -cd ~{in_read_file} | split -l $CHUNK_LINES --filter='pigz -p ~{in_split_read_cores} > ${FILE}.fq.gz' - "fq_chunk_~{in_pair_id}.part."
    
    >>>
    output {
        Array[File] output_read_chunks = glob("fq_chunk_~{in_pair_id}.part.*")
    }
    runtime {
        preemptible: 2
        time: 120
        cpu: in_split_read_cores
        memory: "2 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: "quay.io/glennhickey/pigz:2.3.1"
    }
}

################################################################
# DISTRIBUTE VG-GIRAFFE MAPPING OPERATION OVER EACH CHUNKED READ PAIR #
################################################################

task RUN_VGGIRAFFE {
    input {
        File fastq_file_1
        File? fastq_file_2
        File in_gbz_file
        File in_dist_file
        File in_min_file
        String in_giraffe_options
        String in_sample_name
        Int nb_cores = 16
        String mem_gb = 62
        Int disk_size = 3 * round(size(fastq_file_1, 'G') + size(fastq_file_2, 'G') + size(in_gbz_file, 'G') + size(in_dist_file, 'G') + size(in_min_file, 'G')) + 50
    }
    
    String out_prefix = sub(sub(sub(basename(fastq_file_1), "\\.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")
    Boolean paired_reads = defined(fastq_file_2)
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        PAIR_ARGS=""
        if [ ~{paired_reads} == true ]
        then
            PAIR_ARGS="-f ~{fastq_file_2}"
        fi
        
        vg giraffe \
          --progress \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --sample "~{in_sample_name}" \
          ~{in_giraffe_options} \
          --output-format gaf \
          -f ~{fastq_file_1} ${PAIR_ARGS} \
          -Z ~{in_gbz_file} \
          -d ~{in_dist_file} \
          -m ~{in_min_file} \
          -t ~{nb_cores} | gzip > ~{out_prefix}.gaf.gz
    >>>
    output {
        File chunk_gaf_file = "~{out_prefix}.gaf.gz"
    }
    runtime {
        preemptible: 2
        time: 400
        memory: mem_gb + " GB"
        cpu: nb_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/vgteam/vg:v1.53.0"
    }
}

################################################################
# PREPROCESSING GAF/BAM FILES #
################################################################

task SURJECT_GAF_TO_BAM {
    input {
        File in_gaf_file
        File in_gbz_file
        File in_path_list_file
        String in_sample_name
        Boolean in_paired_reads = true
        Int in_max_fragment_length = 3000
        Boolean input_is_gam = false
        Int thread_count = 16
        String mem_gb = 62
        Int in_disk = round(5 * size(in_gbz_file, 'G') + size(in_gaf_file, 'G')) + 50
    }
    String out_prefix = sub(sub(sub(basename(in_gaf_file), "\\.gz$", ""), "\\.gaf$", ""), "\\.gam$", "")

    command <<<
        set -eux -o pipefail

        PAIR_ARGS=""
        if [ ~{in_paired_reads} == true ]
        then
            PAIR_ARGS="--interleaved --max-frag-len ~{in_max_fragment_length}"
        fi
        
        vg surject \
          -F ~{in_path_list_file} \
          -x ~{in_gbz_file} \
          -t ~{thread_count} \
          --bam-output ~{true="" false="--gaf-input" input_is_gam} \
          --sample ~{in_sample_name} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --prune-low-cplx $PAIR_ARGS \
          ~{in_gaf_file} > ~{out_prefix}.bam
    
    >>>
    output {
        File output_bam_file = "~{out_prefix}.bam"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: mem_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.53.0"
    }
}

task SORT_BAM {
    input {
        File in_bam_file
        File? in_ref_dict
        String in_prefix_to_strip = ""
        Int thread_count = 16
        Int in_disk = round(5 * size(in_bam_file, "G")) + 20
        Int in_mem = 62
    }
    String out_prefix = basename(in_bam_file, ".bam")

    command <<<
        set -eux -o pipefail

        if [ ~{in_prefix_to_strip} != "" ]
        then
            # patch the SQ fields from the dict into a new header
            samtools view -H ~{in_bam_file} | grep ^@HD > new_header.sam
            grep ^@SQ ~{in_ref_dict} | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
            samtools view -H ~{in_bam_file}  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam
            
            cat <(cat new_header.sam) <(samtools view ~{in_bam_file}) | \
                sed -e "s/~{in_prefix_to_strip}//g" | \
                samtools sort --threads ~{thread_count} -O BAM > ~{out_prefix}.positionsorted.bam
        else
            samtools sort --threads ~{thread_count} ~{in_bam_file} \
                     -O BAM > ~{out_prefix}.positionsorted.bam
            
        fi

        samtools index -b ~{out_prefix}.positionsorted.bam ~{out_prefix}.positionsorted.bam.bai
    
    >>>
    output {
        File sorted_bam = "~{out_prefix}.positionsorted.bam"
        File sorted_bam_index = "~{out_prefix}.positionsorted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    }
}

task SPLIT_BAM_BY_PATH {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        String in_prefix_to_strip = ""
        Int thread_count = 16
        Int in_disk = round(3 * size(in_merged_bam_file, "G")) + 20
        Int mem_gb = 62
        Int timeMinutes = 1 + ceil(size(in_merged_bam_file, "G"))
        
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        if [ ~{in_prefix_to_strip} != "" ]
        then
            sed -e "s/~{in_prefix_to_strip}//g" ~{in_path_list_file} > paths.txt
        else
            cp ~{in_path_list_file} paths.txt
        fi
        
        while read -r contig; do
            samtools view \
              -@ ~{thread_count} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < paths.txt

        ## get unmapped reads
        mkdir unmapped
        samtools view \
                 -@ ~{thread_count} \
                 -h -O BAM \
                 -f 4 \
                 input_bam_file.bam \
                 -o unmapped/~{in_sample_name}.unmapped.bam \
    
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
        File bam_unmapped_file = glob("unmapped/~{in_sample_name}.*.bam")[0]
    }
    runtime {
        preemptible: 2
        time: timeMinutes
        memory: mem_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    }
}

task LEFT_SHIFT_BAM_FILE {
    input {
        File in_bam_file
        File in_reference_file
        File in_reference_index_file
        Int in_disk = round(3 * size(in_bam_file, "G")) + 50
        Int in_mem = 62
    }
    String out_prefix = basename(in_bam_file, ".bam")
    
    command <<<
        set -eux -o pipefail 

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        
        bamleftalign \
            < ~{in_bam_file} \
            > ~{out_prefix}.left_shifted.bam \
            --fasta-reference reference.fa \
            --compressed
        samtools index -b ~{out_prefix}.left_shifted.bam ~{out_prefix}.left_shifted.bam.bai
    
    >>>
    output {
        File output_bam_file = "~{out_prefix}.left_shifted.bam"
        File output_bam_index_file = "~{out_prefix}.left_shifted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: in_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/jmonlong/freebayes-samtools:1.2.0_1.10"
    }
}

task PREPARE_REALIGN_TARGETS {
    input {
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_expansion_bases
        Int in_disk = round(2 * size(in_bam_file, "G")) + 20
        Int thread_count = 16
        Int in_mem = 62     
    }
    String out_prefix = basename(in_bam_file, ".bam")
    
    command <<<
        set -eux -o pipefail

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai

        CONTIG_ID=`head -1 < <(samtools view input_bam_file.bam) | cut -f3`
        
        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        # And the dict must be adjacent to both
        ln -f -s "~{in_reference_dict_file}" reference.dict

        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt ~{thread_count} \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{out_prefix}.intervals.bed

        if [ ~{in_expansion_bases} -gt 0 ]; then
            bedtools slop -i ~{out_prefix}.intervals.bed -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > ~{out_prefix}.intervals.widened.bed
            mv ~{out_prefix}.intervals.widened.bed ~{out_prefix}.intervals.bed
        fi
    
    >>>
    output {
        File output_target_bed_file = "~{out_prefix}.intervals.bed"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/jmonlong/gatk-bedtools:3.8.1_2.21.0"
    }
}

task RUN_ABRA_REALIGNER {
    input {
        File in_bam_file
        File in_bam_index_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
        Int in_mem = 62
        Int in_disk = round(3 * (size(in_bam_file, "G") + size(in_reference_file, "G"))) + 50
        Int thread_count = 16
    }
    String out_prefix = basename(in_bam_file, ".bam")
    
    command <<<
        set -eux -o pipefail

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        java '-Xmx~{in_mem}G' -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{out_prefix}.indel_realigned.bam \
          --ref reference.fa \
          --index \
          --threads ~{thread_count}
    
    >>>
    output {
        File indel_realigned_bam = "~{out_prefix}.indel_realigned.bam"
        File indel_realigned_bam_index = "~{out_prefix}.indel_realigned.bai"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task SAMTOOLS_STATS {
    input {
        File in_bam_file        
        Int in_mem = 62       
        Int thread_count = 1
        Int in_disk = round(3 * size(in_bam_file, "G")) + 10
        Int timeMinutes = 1 + ceil(size(in_bam_file, "G"))
    }
    String out_prefix = basename(in_bam_file, ".bam")

    command <<<
        samtools stats --threads ~{thread_count} ~{in_bam_file} > ~{out_prefix}.stats

    >>> 
    output {
        File output_samtools_stats = "~{out_prefix}.stats"
    }
    runtime {
        time: timeMinutes
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"        
        docker: "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    }
}

task MERGE_ALIGNMENT_BAM_CHUNKS {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
        Int in_cores = 16
        Int in_disk = round(5 * size(in_alignment_bam_chunk_files, "G")) + 20
        Int in_mem = 62
    }

    command <<<
        set -eux -o pipefail

        samtools merge \
          -f -p -c --threads ~{in_cores} \
          ~{in_sample_name}_merged.positionsorted.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.positionsorted.bam
    
    >>>
    output {
        File merged_bam_file = "~{in_sample_name}_merged.positionsorted.bam"
        File merged_bam_file_index = "~{in_sample_name}_merged.positionsorted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task MERGE_GAF {
    input {
        String in_sample_name
        Array[File] in_gaf_chunk_files
        Int in_disk = round(3*size(in_gaf_chunk_files, 'G')) + 20
        Int in_mem = 62
    }

    command <<<
        set -eux -o pipefail

        cat ~{sep=" " in_gaf_chunk_files} > ~{in_sample_name}.gaf.gz

    >>>
    output {
        File output_merged_gaf = "~{in_sample_name}.gaf.gz"
    }
    runtime {
        preemptible: 2
        time: 600
        memory: in_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/vgteam/vg:v1.53.0"
    }
}

################################################################
# CALL VARIANTS IN EACH CONTIG #
################################################################

task RUN_DEEP_VARIANT_MAKE_EXAMPLES {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        Int in_min_mapq
        Boolean in_keep_legacy_ac
        Boolean in_norm_reads
        String in_other_makeexamples_arg = ""
        Int in_call_cores
        Int in_call_mem
        Int timeMinutes = 5000
        String in_dv_container = "google/deepvariant:1.5.0"
    }
    Int disk_size = round(2 * size(in_bam_file, 'G')) + 20
    command <<<
        set -eux -o pipefail
        
        ln -s ~{in_bam_file} input_bam_file.bam
        ln -s ~{in_bam_file_index} input_bam_file.bam.bai
        # Files may or may not be indel realigned or left shifted in the names.
        # TODO: move tracking of contig ID to WDL variables!
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
                
        NORM_READS_ARG=""
        if [ ~{in_norm_reads} == true ]; then
          NORM_READS_ARG="--normalize_reads"
        fi

        KEEP_LEGACY_AC_ARG=""
        if [ ~{in_keep_legacy_ac} == true ]; then
          KEEP_LEGACY_AC_ARG="--keep_legacy_allele_counter_behavior"
        fi

        seq 0 $((~{in_call_cores}-1)) | \
        parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref reference.fa \
        --reads input_bam_file.bam \
        --examples ./make_examples.tfrecord@~{in_call_cores}.gz \
        --sample_name ~{in_sample_name} \
        --gvcf ./gvcf.tfrecord@~{in_call_cores}.gz \
        --channels insert_size \
        --min_mapping_quality ~{in_min_mapq} \
        ${KEEP_LEGACY_AC_ARG} ${NORM_READS_ARG} ~{in_other_makeexamples_arg} \
        --regions ${CONTIG_ID} \
        --task {}
        ls | grep 'make_examples.tfrecord-' | tar -czf 'make_examples.tfrecord.tar.gz' -T -
        ls | grep 'gvcf.tfrecord-' | tar -czf 'gvcf.tfrecord.tar.gz' -T -
    >>>
    output {
        File examples_file = "make_examples.tfrecord.tar.gz"
        File nonvariant_site_tf_file = "gvcf.tfrecord.tar.gz"
    }
    runtime {
        preemptible: 5
        time: timeMinutes
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        disks: "local-disk " + disk_size + " SSD"
        docker: in_dv_container
    }
}

task RUN_DEEP_VARIANT {
    input {
        String in_sample_name
        File in_reference_file
        File in_reference_index_file
        File in_examples_file
        File in_nonvariant_site_tf_file
        File? in_model_meta_file
        File? in_model_index_file
        File? in_model_data_file
        Boolean? VCFStatsReport = true
        Int in_call_cores
        Int in_call_mem
        Int timeMinutes = 5000
        String in_dv_gpu_container = "google/deepvariant:1.5.0-gpu"
    }
    Int disk_size = 5 * round(size(in_examples_file, 'G') + size(in_nonvariant_site_tf_file, 'G') + size(in_reference_file, 'G')) + 50
    command <<<
        set -eux -o pipefail
        
        tar -xzf ~{in_examples_file}
        tar -xzf ~{in_nonvariant_site_tf_file}
        
        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        # We should use an array here, but that doesn't seem to work the way I
        # usually do them (because of a set -u maybe?)
        if [[ ! -z "~{in_model_meta_file}" ]] ; then
            # Model files must be adjacent and not at arbitrary paths
            ln -f -s "~{in_model_meta_file}" model.meta
            ln -f -s "~{in_model_index_file}" model.index
            ln -f -s "~{in_model_data_file}" model.data-00000-of-00001
        else
            # use default models
            ln -f -s "/opt/models/wes/model.ckpt.meta" model.meta
            ln -f -s "/opt/models/wes/model.ckpt.index" model.index
            ln -f -s "/opt/models/wes/model.ckpt.data-00000-of-00001" model.data-00000-of-00001
        fi
        
        /opt/deepvariant/bin/call_variants \
        --outfile call_variants_output.tfrecord.gz \
        --examples "make_examples.tfrecord@~{in_call_cores}.gz" \
        --checkpoint model && \
        /opt/deepvariant/bin/postprocess_variants \
        --ref reference.fa \
        --infile call_variants_output.tfrecord.gz \
        --nonvariant_site_tfrecord_path "gvcf.tfrecord@~{in_call_cores}.gz" \
        --outfile "~{in_sample_name}_deepvariant.vcf.gz" \
        --gvcf_outfile "~{in_sample_name}_deepvariant.g.vcf.gz" \
        ~{true="--vcf_stats_report" false="--novcf_stats_report" VCFStatsReport} \
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deepvariant.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deepvariant.g.vcf.gz"
        Array[File] outputVCFStatsReport = glob("*.visual_report.html")
    }
    runtime {
        preemptible: 5
        time: timeMinutes
        maxRetries: 5
        memory: in_call_mem + " GB"
        cpu: in_call_cores
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + disk_size + " SSD"
        docker: in_dv_gpu_container
    }
}

task CONCAT_CLIPPED_VCF_CHUNKS {
    input {
        String in_sample_name
        String in_tools
        Array[File] in_clipped_vcf_chunk_files
        Int in_mem = 8
        Int in_disk = round(30 * size(in_clipped_vcf_chunk_files, "G")) + 50      
    }

    command {
        set -eux -o pipefail

        mkdir bcftools.tmp
        bcftools concat -n ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort -T bcftools.tmp -O z -o ${in_sample_name}.${in_tools}.vcf.gz - && bcftools index -t -o ${in_sample_name}.${in_tools}.vcf.gz.tbi ${in_sample_name}.${in_tools}.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}.${in_tools}.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}.${in_tools}.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        time: 60
        memory: in_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    }
}

################################################################
# LIFT GENOME ANNOTATIONS #
################################################################

task EXTRACT_TARGET_REGION_BAM {
    input {
        File in_bam_file 
        File in_bam_index_file
        File in_new_header_file
        String bam_target_region      
        Int thread_count = 1
        Int in_mem = 4
        Int timeMinutes = 1 + ceil(size(in_bam_file, "G"))
        Int in_disk = round(5 * size(in_bam_file, "G")) + 20
    }
    String out_prefix = basename(in_bam_file, ".bam")

    command <<<

        samtools view -h -o ~{out_prefix}.target.bam ~{in_bam_file} ~{bam_target_region}
        
        samtools index ~{out_prefix}.target.bam

        cp ~{in_new_header_file} new_header.tsv
        
        samtools view -H ~{out_prefix}.target.bam | grep -e "^@RG" >> new_header.tsv        
              
        samtools reheader new_header.tsv ~{out_prefix}.target.bam > ~{out_prefix}.target.reheader.bam
        
        samtools index ~{out_prefix}.target.reheader.bam

    >>> 
    output {
        File target_bam = "~{out_prefix}.target.reheader.bam"
        File target_bam_index = "~{out_prefix}.target.reheader.bam.bai"
    }
    runtime {
        time: timeMinutes
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"        
        docker: "quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    }
}

task LIFTOVER_BAM_FILE {
    input {
        File in_bam_file 
        File in_bam_index_file
        File in_position_report
        File script               
        Int in_mem = 62
        Int in_disk = round(5 * size(in_bam_file, "G")) + 50
        Int thread_count = 16
    }
    String out_prefix = basename(in_bam_file, ".bam")

    command <<<
        set -eux -o pipefail          
        python ~{script} ~{in_position_report} ~{in_bam_file} ~{out_prefix}.liftover.bam  

    >>>
    output {
        File liftover_bam = "~{out_prefix}.liftover.sorted.bam"
        File liftover_bam_index = "~{out_prefix}.liftover.sorted.bam.bai"        
    }
    runtime {
        preemptible: 2
        time: 9000
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/pysam:0.22.0--py38h15b938a_0"
    }
}

task LIFTOVER_VCF_FILE {
    input {
        File in_vcf_file
        File in_position_report
        File script               
        Int in_mem = 16
        Int in_disk = round(2 * size(in_vcf_file, "G")) + 10
        Int timeMinutes = 1 + ceil(size(in_vcf_file, "G"))
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")

    command <<<
        set -eux -o pipefail       
        python ~{script} ~{in_position_report} ~{in_vcf_file} ~{out_prefix}.liftover.vcf  

    >>>
    output {
        File vcf_liftover = "~{out_prefix}.liftover.vcf"
    }
    runtime {
        preemptible: 2
        time: timeMinutes
        cpu: 8
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/python:3.10.2"
    }
}

task REHEADER_VCF {
    input {
        File in_vcf_file
        File in_reference_index_file
        Int in_mem = 1
        Int in_disk = round(2 * size(in_vcf_file, "G")) + 10
        Int timeMinutes = 1 + ceil(size(in_vcf_file, "G"))     
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")

    command <<<
        set -eux -o pipefail

        mkdir bcftools.tmp
        
        bcftools reheader -f ~{in_reference_index_file} ~{in_vcf_file} -o ~{out_prefix}.original.reheader.vcf
        
        bcftools sort -T bcftools.tmp -O z -o ~{out_prefix}.reheader.vcf.gz ~{out_prefix}.original.reheader.vcf - && bcftools index -t -o ~{out_prefix}.reheader.vcf.gz.tbi ~{out_prefix}.reheader.vcf.gz
        
    >>>
    output {
        File reheader_vcf = "~{out_prefix}.reheader.vcf.gz"
        File reheader_vcf_index = "~{out_prefix}.reheader.vcf.gz.tbi"       
    }
    runtime {
        preemptible: 2
        time: timeMinutes
        cpu: 1
        memory: in_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    } 
}

################################################################
# FILTER VARIANTS #
################################################################

task VCF_FILTER {
    input {
        File in_vcf_file
        File in_vcf_file_index
        Int in_mem = 16
        Int in_disk = round(3 * size(in_vcf_file, "G")) + 20
        Int timeMinutes = 1 + ceil(size(in_vcf_file, "G")) * 5
        Int in_cores = 8   
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")

    command <<<
        set -eux -o pipefail

        rtg vcffilter -i ~{in_vcf_file} -o ~{out_prefix}.filtered.vcf.gz --min-quality=20 --min-read-depth=10
    
    >>>
    output {
        File filtered_vcf = "~{out_prefix}.filtered.vcf.gz"
        File filtered_vcf_index = "~{out_prefix}.filtered.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        time: timeMinutes
        cpu: in_cores
        memory: in_mem + " GB"             
        disks: "local-disk " + in_disk + " SSD"
        docker: "realtimegenomics/rtg-tools:3.12.1"       
    }  
}

################################################################
# VCF STATISTIC #
################################################################

task BCFTOOLS_STATS {
    input {
        File in_vcf_file
        File in_vcf_file_index
        Int in_mem = 4        
        Int thread_count = 1
        Int in_disk = round(2 * size(in_vcf_file, "G")) + 10
        Int timeMinutes = 1 + ceil(size(in_vcf_file, "G")) * 2   
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")

    command <<<
        bcftools stats --threads ~{thread_count} ~{in_vcf_file} > ~{out_prefix}.bcftools_stats.txt

    >>> 
    output {
        File bcftools_stats = "~{out_prefix}.bcftools_stats.txt"
    }
    runtime {
        time: timeMinutes
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"        
        docker: "quay.io/biocontainers/bcftools:1.19--h8b25389_0"
    }
}

task VCFTOOLS_SUMMARY {
    input {
        File in_vcf_file
        File in_vcf_file_index
        Int in_mem = 4        
        Int thread_count = 1
        Int in_disk = round(2 * size(in_vcf_file, "G")) + 10
        Int timeMinutes = 1 + ceil(size(in_vcf_file, "G")) * 2  
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")
    
    command <<<
        vcftools --TsTv-by-count --gzvcf ~{in_vcf_file} --out ~{out_prefix}
        vcftools --TsTv-by-qual --gzvcf ~{in_vcf_file} --out ~{out_prefix}
        vcftools --FILTER-summary --gzvcf ~{in_vcf_file} --out ~{out_prefix}

    >>> 
    output {
        File TsTv_count = "~{out_prefix}.TsTv.count"
        File TsTv_qual = "~{out_prefix}.TsTv.qual"
        File FILTER_summary = "~{out_prefix}.FILTER.summary"
    }
    runtime {
        time: timeMinutes
        memory: in_mem + " GB"
        cpu: thread_count
        disks: "local-disk " + in_disk + " SSD"        
        docker: "quay.io/biocontainers/vcftools:0.1.16--he513fc3_4"
    }
}

################################################################
# SNPEFF: ANNOTATIONS #
################################################################

task SNPEFF_ANNOTATE_VCF {
    input {
        File in_vcf_file
        File in_vcf_file_index
        File? in_snpeff_database
        Int in_call_cores
        Int in_call_mem
        Int in_disk = round(3 * (size(in_vcf_file, "G") + size(in_snpeff_database, "G"))) + 50
        Int timeMinutes = 400    
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")
    
    command <<<
        set -eux -o pipefail
        
        unzip ~{in_snpeff_database}
        snpEff '-Xmx~{in_call_mem}G' -i VCF -o VCF -noLof -noHgvs -formatEff -classic -dataDir ${PWD}/data GRCh38.105 ~{in_vcf_file} -csvStats ~{out_prefix}.stats.csv > ~{out_prefix}.snpeff.vcf
    
    >>>
    output {
        File annotated_vcf = "~{out_prefix}.snpeff.vcf"
        File annotated_stats = "~{out_prefix}.stats.csv"
        File annotated_html = "snpEff_summary.html"
        File annotated_gene = "~{out_prefix}.stats.genes.txt"
    }
    runtime {
       time: timeMinutes
        cpu: in_call_cores
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_disk + " SSD"
        docker: "quay.io/biocontainers/snpeff:5.1--hdfd78af_2"
    }
}

task SNPEFF_ANNOTATE_VCF_LOCAL {
    input {
        File in_vcf_file
        File in_vcf_file_index
        String? in_name_snpeff_database
        Int in_mem = 62
        Int in_cores = 16
        Int in_disk = round(3 * (size(in_vcf_file, "G"))) + 50
        Int timeMinutes = 120
    }
    String out_prefix = basename(in_vcf_file, ".vcf.gz")

    command <<<
        set -euo pipefail

        # Set the path to the snpEff.jar and SnpSift.jar files on your local machine
        snpEffJarPath="$HOME/snpEff/snpEff.jar"
        snpSiftJarPath="$HOME/snpEff/SnpSift.jar"
        
        # Define the path to the ClinVar VCF file
        clinvar="$HOME/snpEff/~{in_name_snpeff_database}/~{in_name_snpeff_database}.clinvar.vcf.gz"
        
        # Run snpEff to annotate genetic variants
        java '-Xmx~{in_mem}G' -jar $snpEffJarPath ~{in_name_snpeff_database} ~{in_vcf_file} -csvStats ~{out_prefix}.stats.csv > ~{out_prefix}.ann.vcf

        gatk --java-options '-Xmx~{in_mem}G' VariantsToTable -V ~{out_prefix}.ann.vcf -F CHROM -F POS -F TYPE -F ID -F ANN -F LOF -F NMD -GF AD -GF DP -GF GQ -GF GT -O ~{out_prefix}.ann.vcf.csv

        bgzip ~{out_prefix}.ann.vcf && tabix -p vcf ~{out_prefix}.ann.vcf.gz

        mv snpEff_summary.html ~{out_prefix}.snpEff_summary.html

        # Run snpEff to annotate genetic variants with ClinVar
        java  '-Xmx~{in_mem}G' -jar $snpSiftJarPath annotate $clinvar ~{in_vcf_file} > ~{out_prefix}.ann.clinvar.vcf

        gatk --java-options '-Xmx~{in_mem}G' VariantsToTable -V ~{out_prefix}.ann.clinvar.vcf -F CHROM -F POS -F TYPE -F ID -F ALLELEID -F CLNDN -F CLNSIG -F CLNSIGCONF -F CLNSIGINCL -F CLNVC -F GENEINFO -GF AD -GF GQ -GF GT -O ~{out_prefix}.ann.clinvar.vcf.csv 

        bgzip ~{out_prefix}.ann.clinvar.vcf && tabix -p vcf ~{out_prefix}.ann.clinvar.vcf.gz

        # Extract variants annotated as either Pathogenic or Likely_pathogenic.
        awk -F'\t' '$7 ~ /^Pathogenic|^Likely_pathogenic/' ~{out_prefix}.ann.clinvar.vcf.csv > ~{out_prefix}.ann.clinvar.pathogenic_list.csv

        if [ -s ~{out_prefix}.ann.clinvar.pathogenic_list.csv ]; then
            awk -F'\t' '{print $1 "\t" ($2 - 1) "\t" $2}' ~{out_prefix}.ann.clinvar.pathogenic_list.csv > ~{out_prefix}.ann.clinvar.pathogenic_list.interval.bed    
            rtg vcffilter -i ~{in_vcf_file} -o ~{out_prefix}.ann.clinvar.pathogenic.vcf.gz --include-bed ~{out_prefix}.ann.clinvar.pathogenic_list.interval.bed
        fi    
        
        if [ ! -s ~{out_prefix}.ann.clinvar.pathogenic_list.csv ]; then
            echo "There are no variants identified as pathogenic."
        fi

    >>>
    output {
        File annotated_vcf = "~{out_prefix}.ann.vcf.gz"
        File annotated_vcf_index = "~{out_prefix}.ann.vcf.gz.tbi"
        File annotated_stats = "~{out_prefix}.stats.csv"
        File annotated_csv_report = "~{out_prefix}.ann.vcf.csv"
        File annotated_clinvar_vcf = "~{out_prefix}.ann.clinvar.vcf.gz"
        File annotated_clinvar_vcf_index = "~{out_prefix}.ann.clinvar.vcf.gz.tbi"
        File annotated_clinvar_csv_report = "~{out_prefix}.ann.clinvar.vcf.csv"
        File? annotated_pathogenic_list = "~{out_prefix}.ann.clinvar.pathogenic_list.csv"
        File? annotated_pathogenic_interval = "~{out_prefix}.ann.clinvar.pathogenic_list.interval.bed"
        File? annotated_pathogenic_vcf = "~{out_prefix}.ann.clinvar.pathogenic.vcf.gz"
        File? annotated_pathogenic_vcf_index = "~{out_prefix}.ann.clinvar.pathogenic.vcf.gz.tbi"
        File annotated_html = "~{out_prefix}.snpEff_summary.html"
        File annotated_gene = "~{out_prefix}.stats.genes.txt"
    }
}

################################################################
# MULTIQC: SUMMARIZE BIOINFORMATICS FINDINGS FROM MULTIPLE SAMPLES INTO ONE REPORT #
################################################################

task MULTIQC {
    input {
        Array[File] input_files
        Array[File] samtools_stats
        String prefix
        File config
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float input_size = size(input_files, "G")
    Int disk_size_gb = ceil(input_size) + 5 + modify_disk_size_gb

    String out_tar_gz = prefix + ".tar.gz"

    command <<<
        set -euo pipefail

        cp -r ../inputs .
        #cp -r _miniwdl_inputs/0 .
        
        multiqc -v --force --config ~{config} -o ~{prefix} .

        if [ ! -d ~{prefix} ]; then
            >&2 echo "MultiQC didn't find any valid files!"
            exit 1
        fi

        tar -czf ~{out_tar_gz} ~{prefix}

    >>>
    output {
        File multiqc_report = out_tar_gz
    }
    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'
        maxRetries: max_retries
    }
}
        
