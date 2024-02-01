version 1.0

import "path/to/biotask_utils.wdl" as utils

workflow UP_DOWNSTREAM_PANGENOME {

    meta {
        description: "## UP/DOWNSTREAM_PANGENOME workflow. \n The full workflow to go from sequencing reads (FASTQs, CRAM) to small variant calls (VCF). Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. left-align, indel realignment). DeepVariant calls germline small variants. Finally, SnpEff annotates the identified variants. \n This WDL workflow draws its foundation from the work of Liao, WW., Asri, M., Ebler, J. et al. A draft human pangenome reference. Nature 617, 312–324 (2023). [https://doi.org/10.1038/s41586-023-05896-x] (https://github.com/vgteam/vg_wdl)."
    }   

    parameter_meta {
        INPUT_READ_FILE_1: "Input sample 1st read pair fastq.gz"
        INPUT_READ_FILE_2: "Input sample 2nd read pair fastq.gz"
        INPUT_CRAM_FILE: "Input CRAM file"
        CRAM_REF: "Genome fasta file associated with the CRAM file"
        CRAM_REF_INDEX: "Index of the fasta file associated with the CRAM file"
        GBZ_FILE: "Path to .gbz index file"
        DIST_FILE: "Path to .dist index file"
        MIN_FILE: "Path to .min index file"
        SAMPLE_NAME: "The sample name"
        OUTPUT_GAF: "Should a GAF file with the aligned reads be saved? Default is 'true'."
        OUTPUT_SINGLE_BAM: "Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling bams' (one per contig) won't be outputed. Default is 'true'."
        PAIRED_READS: "Are the reads paired? Default is 'true'."
        TRIM_FASTQ: "Run FastP for read trimming? Default is 'true'."
        READS_PER_CHUNK: "Number of reads contained in each mapping chunk. Default 20 000 000."
        PATH_LIST_FILE: "(OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths."
        CONTIGS: "(OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index."
        REFERENCE_PREFIX: "Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)"
        REFERENCE_FILE: "(OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference."
        REFERENCE_INDEX_FILE: "(OPTIONAL) If specified, use this .fai index instead of indexing the reference file."
        REFERENCE_DICT_FILE: "(OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set"
        LEFTALIGN_BAM: "Whether or not to left-align reads in the BAM. Default is 'true'."
        REALIGN_INDELS: "Whether or not to realign reads near indels. Default is 'true'."
        REALIGNMENT_EXPANSION_BASES: "Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160."
        MIN_MAPQ: "Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is 1"
        MAX_FRAGMENT_LENGTH: "Maximum distance at which to mark paired reads properly paired. Default is 3000."
        GIRAFFE_OPTIONS: "(OPTIONAL) extra command line options for Giraffe mapper"
        DV_MODEL_META: ".meta file for a custom DeepVariant calling model"
        DV_MODEL_INDEX: ".index file for a custom DeepVariant calling model"
        DV_MODEL_DATA: ".data-00000-of-00001 file for a custom DeepVariant calling model"
        DV_KEEP_LEGACY_AC: "Should DV use the legacy allele counter behavior? Default is 'true'."
        DV_NORM_READS: "Should DV normalize reads itself? Default is 'fasle'."
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant."
	LIFTOVER_BAM: "When dealing with a pangenome that has different chromosome names compared to the standard linear reference genome, converts genomic BAM data between reference assemblies. Default is 'false'."
        LIFTOVER_BAM_SCRIPT: "A script converts VCF Ffiles between different genome builds. Required if using LIFTOVER_BAM."
        LIFTOVER_BAM_TARGET_REGION: "Liftover BAM file is solely for IGV visualization, and this process is time-consuming. Therefore, you should extract the specified target regions.Required if using LIFTOVER_BAM."
        LIFTOVER_BAM_NEW_HEADER: "Pangenome BAM header differs, submit linear pipeline's header. Required if using LIFTOVER_BAM."
        LIFTOVER_VCF: "When dealing with a pangenome that has different chromosome names compared to the standard linear reference genome, converts genomic VCF data between reference assemblies. Default is 'false'."
        LIFTOVER_VCF_SCRIPT: "A script converts VCF Ffiles between different genome builds. Required if using LIFTOVER_VCF."
        LIFTOVER_POSITION_REPORT: "A TSV file that records the ordinal positions of the standard linear reference genome. Required if using LIFTOVER_VCF." 
        LIFTOVER_LINEAR_REFERENCE_INDEX: "Utilize the .fai index of the standard linear reference genome to reconfigure lifted-over VCF files. Required if using LIFTOVER_VCF."
        SNPEFF_ANNOTATION: "Set to 'true' to run SnpEff annotation on the joint genotyped VCF on docker."
        SNPEFF_DATABASE: "Path to SnpEff database .zip file for SnpEff annotation functionality. Required if using SNPEFF_ANNOTATION."
        SNPEFF_ANNOTATION_LOCAL: "Set to 'true' to run SnpEff and ClinVar snpSift annotation on your local machine for faster processing. Be sure to specify the paths for the SnpEff.jar, SnpSift.jar, and ClinVar VCF file. Default is 'false'"
        NAME_SNPEFF_DATABASE: "Name of SnpEff genome version. Required if using SNPEFF_ANNOTATION_LOCAL."
        MULTIQC_CONFIG: "YAML file for configuring MultiQC"
        SPLIT_READ_CORES: "Number of cores to use when splitting the reads into chunks. Default is 8."
        MAP_CORES: "Number of cores to use when mapping the reads. Default is 16."
        MAP_MEM: "Memory, in GB, to use when mapping the reads. Default is 120."
        CALL_CORES: "Number of cores to use when calling variants. Default is 8."
        CALL_MEM: "Memory, in GB, to use when calling variants. Default is 50."
    }

    input {
        File? INPUT_READ_FILE_1
        File? INPUT_READ_FILE_2
        File? INPUT_CRAM_FILE
        File? CRAM_REF
        File? CRAM_REF_INDEX
        File GBZ_FILE
        File DIST_FILE
        File MIN_FILE
        String SAMPLE_NAME
        Boolean OUTPUT_GAF = true
        Boolean OUTPUT_SINGLE_BAM = true
        Boolean PAIRED_READS = true
        Boolean TRIM_FASTQ = true
        Int READS_PER_CHUNK = 20000000
        File? PATH_LIST_FILE
        Array[String]+? CONTIGS
        String REFERENCE_PREFIX = ""
        File? REFERENCE_FILE
        File? REFERENCE_INDEX_FILE
        File? REFERENCE_DICT_FILE                
        Boolean LEFTALIGN_BAM = true
        Boolean REALIGN_INDELS = true
        Int REALIGNMENT_EXPANSION_BASES = 160
        Int MIN_MAPQ = 1
        Int MAX_FRAGMENT_LENGTH = 3000
        String GIRAFFE_OPTIONS = ""
        File? DV_MODEL_META
        File? DV_MODEL_INDEX
        File? DV_MODEL_DATA
        Boolean DV_KEEP_LEGACY_AC = true
        Boolean DV_NORM_READS = false
        String OTHER_MAKEEXAMPLES_ARG = ""
        Boolean LIFTOVER_BAM = false
        File? LIFTOVER_BAM_SCRIPT
        String LIFTOVER_BAM_TARGET_REGION = ""
        File? LIFTOVER_BAM_NEW_HEADER
        Boolean LIFTOVER_VCF = false
        File? LIFTOVER_VCF_SCRIPT
        File? LIFTOVER_POSITION_REPORT
        File? LIFTOVER_LINEAR_REFERENCE_INDEX
        Boolean SNPEFF_ANNOTATION = false
        File? SNPEFF_DATABASE
        Boolean SNPEFF_ANNOTATION_LOCAL = false
        String? NAME_SNPEFF_DATABASE = ""
        File? MULTIQC_CONFIG
        Int SPLIT_READ_CORES = 8
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Int CALL_CORES = 16
        Int CALL_MEM = 50
    }

################################################################
# PREPARE DATA #
################################################################

    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from GBZ file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call utils.EXTRACTSUBSETPATHNAMES {
                input:
                    in_gbz_file=GBZ_FILE,
                    in_extract_mem=MAP_MEM
            }
        }
    }

    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }

    File pipeline_path_list_file = select_first([PATH_LIST_FILE, EXTRACTSUBSETPATHNAMES.output_path_list_file, written_path_names_file])
    
    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph, we generate it ourselves, from the graph.
    if (!defined(REFERENCE_FILE)) {
        call utils.EXTRACTREFERENCE {
            input:
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=MAP_MEM
        }
    }

    File reference_file = select_first([REFERENCE_FILE, EXTRACTREFERENCE.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.INDEXREFERENCE {
            input:
                in_reference_file=reference_file
        }
    }

    File reference_index_file = select_first([REFERENCE_INDEX_FILE, INDEXREFERENCE.reference_index_file])
    
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, INDEXREFERENCE.reference_dict_file])
    
    if (defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
	    call utils.CONVERT_CRAM_TO_FASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES
	    }
    }

################################################################
# SEQUENCING QUALITY CONTROL #
################################################################
    
    if (!PAIRED_READS) {
        call utils.FASTQC_SINGLE_END {
            input:
                in_sample_name = SAMPLE_NAME,
                in_read_1_fastq = select_first([INPUT_READ_FILE_1, CONVERT_CRAM_TO_FASTQ.output_fastq_1_file])
        }
        if (TRIM_FASTQ) {
            call utils.FASTP_SINGLE_END {
                input:
                    in_sample_name = SAMPLE_NAME,
                    in_read_1_fastq = select_first([INPUT_READ_FILE_1, CONVERT_CRAM_TO_FASTQ.output_fastq_1_file]),
                    in_call_mem = CALL_MEM
            }
        }        
    } 

    if (PAIRED_READS) {
        call utils.FASTQC_PAIRED_READS {
            input:
                in_sample_name = SAMPLE_NAME,
                in_read_1_fastq = select_first([INPUT_READ_FILE_1, CONVERT_CRAM_TO_FASTQ.output_fastq_1_file]),
                in_read_2_fastq = select_first([INPUT_READ_FILE_2, CONVERT_CRAM_TO_FASTQ.output_fastq_2_file])
        }
        if  (TRIM_FASTQ) {
            call utils.FASTP_PAIRED_READS {
                input:
                    in_sample_name = SAMPLE_NAME,
                    in_read_1_fastq = select_first([INPUT_READ_FILE_1, CONVERT_CRAM_TO_FASTQ.output_fastq_1_file]),
                    in_read_2_fastq = select_first([INPUT_READ_FILE_2, CONVERT_CRAM_TO_FASTQ.output_fastq_2_file]),
                    in_call_mem = CALL_MEM
            }
        }        
    }   
     
################################################################
# DISTRIBUTE VG-GIRAFFE MAPPING OPERATION OVER EACH CHUNKED READ PAIR #
################################################################
    File read_1_file = select_first([FASTP_PAIRED_READS.output_fastp_1_file, FASTP_SINGLE_END.output_fastp_1_file, INPUT_READ_FILE_1, CONVERT_CRAM_TO_FASTQ.output_fastq_1_file])
    
    # Split input reads into chunks for parallelized mapping
    call utils.SPLIT_READS as firstReadPair {
        input:
            in_read_file=read_1_file,
            in_pair_id="1",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
    }

    if (PAIRED_READS) {
        File read_2_file = select_first([FASTP_PAIRED_READS.output_fastp_2_file, INPUT_READ_FILE_2, CONVERT_CRAM_TO_FASTQ.output_fastq_2_file])
        call utils.SPLIT_READS as secondReadPair {
            input:
            in_read_file=read_2_file,
            in_pair_id="2",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
        }
        Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
        scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
            call utils.RUN_VGGIRAFFE as RUN_VGGIRAFFE_PE {
                input:
                fastq_file_1=read_pair_chunk_files.left,
                fastq_file_2=read_pair_chunk_files.right,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=GBZ_FILE,
                in_dist_file=DIST_FILE,
                in_min_file=MIN_FILE,
                # We always need to pass a full dict file here, with lengths,
                # because if we pass just path lists and the paths are not
                # completely contained in the graph (like if we're working on
                # GRCh38 paths in a CHM13-based graph), giraffe won't be able
                # to get the path lengths and will crash.
                # TODO: Somehow this problem is supposed to go away if we pull
                # any GRCh38. prefix off the path names by setting
                # REFERENCE_PREFIX and making sure the prefix isn't in the
                # truth set.
                # See <https://github.com/adamnovak/giraffe-dv-wdl/pull/2#issuecomment-955096920>
                in_sample_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                mem_gb=MAP_MEM
            }
        }
    }

    if (!PAIRED_READS) {
        scatter (read_pair_chunk_file in firstReadPair.output_read_chunks) {
            call utils.RUN_VGGIRAFFE as RUN_VGGIRAFFE_SE {
                input:
                fastq_file_1=read_pair_chunk_file,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_gbz_file=GBZ_FILE,
                in_dist_file=DIST_FILE,
                in_min_file=MIN_FILE,
                # We always need to pass a full dict file here, with lengths,
                # because if we pass just path lists and the paths are not
                # completely contained in the graph (like if we're working on
                # GRCh38 paths in a CHM13-based graph), giraffe won't be able
                # to get the path lengths and will crash.
                # TODO: Somehow this problem is supposed to go away if we pull
                # any GRCh38. prefix off the path names by setting
                # REFERENCE_PREFIX and making sure the prefix isn't in the
                # truth set.
                # See <https://github.com/adamnovak/giraffe-dv-wdl/pull/2#issuecomment-955096920>
                in_sample_name=SAMPLE_NAME,
                nb_cores=MAP_CORES,
                mem_gb=MAP_MEM
            }
        }
    }

    Array[File] gaf_chunks = select_first([RUN_VGGIRAFFE_PE.chunk_gaf_file, RUN_VGGIRAFFE_SE.chunk_gaf_file])
    scatter (gaf_file in gaf_chunks) {
        call utils.SURJECT_GAF_TO_BAM {
            input:
            in_gaf_file=gaf_file,
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_sample_name=SAMPLE_NAME,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            in_paired_reads=PAIRED_READS,
            mem_gb=MAP_MEM
        }

        call utils.SORT_BAM {
            input:
            in_bam_file=SURJECT_GAF_TO_BAM.output_bam_file,
            in_ref_dict=reference_dict_file,
            in_prefix_to_strip=REFERENCE_PREFIX
        }
    }

    call utils.MERGE_ALIGNMENT_BAM_CHUNKS {
        input:
        in_sample_name=SAMPLE_NAME,
        in_alignment_bam_chunk_files=SORT_BAM.sorted_bam
    }

    # Split merged alignment by contigs list
    call utils.SPLIT_BAM_BY_PATH {
        input:
        in_sample_name=SAMPLE_NAME,
        in_merged_bam_file=MERGE_ALIGNMENT_BAM_CHUNKS.merged_bam_file,
        in_merged_bam_file_index=MERGE_ALIGNMENT_BAM_CHUNKS.merged_bam_file_index,
        in_path_list_file=pipeline_path_list_file,
        in_prefix_to_strip=REFERENCE_PREFIX
    }

################################################################
# CALL VARIANTS IN EACH CONTIG #
################################################################  
    
    scatter (bam_and_index_for_path in zip(SPLIT_BAM_BY_PATH.bam_contig_files, SPLIT_BAM_BY_PATH.bam_contig_files_index)) {
        ## Evantually shift and realign reads
        if (LEFTALIGN_BAM) {
            # Just left-shift each read individually
            call utils.LEFT_SHIFT_BAM_FILE {
                input:
                in_bam_file=bam_and_index_for_path.left,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file
            }
        }

        if (REALIGN_INDELS) {
            File forrealign_bam = select_first([LEFT_SHIFT_BAM_FILE.output_bam_file, bam_and_index_for_path.left])
            File forrealign_index = select_first([LEFT_SHIFT_BAM_FILE.output_bam_index_file, bam_and_index_for_path.right])
            # Do indel realignment
            call utils.PREPARE_REALIGN_TARGETS {
                input:
                in_bam_file=forrealign_bam,
                in_bam_index_file=forrealign_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_reference_dict_file=reference_dict_file,
                in_expansion_bases=REALIGNMENT_EXPANSION_BASES
            }
            call utils.RUN_ABRA_REALIGNER {
                input:
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_target_bed_file=PREPARE_REALIGN_TARGETS.output_target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,                    
            }
        }

        File calling_bam = select_first([RUN_ABRA_REALIGNER.indel_realigned_bam, LEFT_SHIFT_BAM_FILE.output_bam_file, bam_and_index_for_path.left])
        File calling_bam_index = select_first([RUN_ABRA_REALIGNER.indel_realigned_bam_index, LEFT_SHIFT_BAM_FILE.output_bam_index_file, bam_and_index_for_path.right])

        call utils.SAMTOOLS_STATS as SAMTOOLS_STATS_CALLING_BAM {
            input:                 
            in_bam_file=calling_bam
        } 

        if (LIFTOVER_BAM) {
            call utils.EXTRACT_TARGET_REGION_BAM {
                input:
                    in_bam_file=calling_bam,
                    in_bam_index_file=calling_bam_index,
                    bam_target_region=LIFTOVER_BAM_TARGET_REGION,
                    in_new_header_file=select_first([LIFTOVER_BAM_NEW_HEADER,[]])
            }

            call utils.LIFTOVER_BAM_FILE {
                input:
                    in_bam_file=EXTRACT_TARGET_REGION_BAM.target_bam,
                    in_bam_index_file=EXTRACT_TARGET_REGION_BAM.target_bam_index,
                    in_position_report=select_first([LIFTOVER_POSITION_REPORT,[]]),
                    script=select_first([LIFTOVER_BAM_SCRIPT,[]])                        
            }           
        }                            
        
        File calling_bam_liftover =select_first([LIFTOVER_BAM_FILE.liftover_bam,bam_and_index_for_path.left])
        File calling_bam_liftover_index =select_first([LIFTOVER_BAM_FILE.liftover_bam_index,bam_and_index_for_path.right])                          
        
        ## DEEPVARIANT calling
        call utils.RUN_DEEP_VARIANT_MAKE_EXAMPLES {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=calling_bam,
                in_bam_file_index=calling_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_min_mapq=MIN_MAPQ,
                in_keep_legacy_ac=DV_KEEP_LEGACY_AC,
                in_norm_reads=DV_NORM_READS,
                in_other_makeexamples_arg=OTHER_MAKEEXAMPLES_ARG,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
        }

        call utils.RUN_DEEP_VARIANT {
            input:
                in_sample_name=SAMPLE_NAME,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_examples_file=RUN_DEEP_VARIANT_MAKE_EXAMPLES.examples_file,
                in_nonvariant_site_tf_file=RUN_DEEP_VARIANT_MAKE_EXAMPLES.nonvariant_site_tf_file,
                in_model_meta_file=DV_MODEL_META,
                in_model_index_file=DV_MODEL_INDEX,
                in_model_data_file=DV_MODEL_DATA,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
        }              
    }

    # Merge distributed variant called VCFs
    call utils.CONCAT_CLIPPED_VCF_CHUNKS as CONCATCLIPPEDVCFCHUNKS_DEEPVARIANT {
        input:
            in_sample_name=SAMPLE_NAME,
            in_tools="deepvariant",
            in_clipped_vcf_chunk_files=RUN_DEEP_VARIANT.output_vcf_file
    }
    call utils.CONCAT_CLIPPED_VCF_CHUNKS as CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT {
        input:
            in_sample_name=SAMPLE_NAME,
            in_tools="deepvariant.g",
            in_clipped_vcf_chunk_files=RUN_DEEP_VARIANT.output_gvcf_file
    }

################################################################
# OUTPUT GAF/BAM FILES #
################################################################     
    
    if (LIFTOVER_BAM) {
        Array[File] output_calling_bam_liftover_files = calling_bam_liftover
        Array[File] output_calling_bam_liftover_index_files = calling_bam_liftover_index  
    }

    if (OUTPUT_SINGLE_BAM) {
        call utils.MERGE_ALIGNMENT_BAM_CHUNKS as MERGE_BAM {
            input:
                in_sample_name=SAMPLE_NAME,
                in_alignment_bam_chunk_files=select_all(flatten([calling_bam, [SPLIT_BAM_BY_PATH.bam_unmapped_file]]))
        }
        call utils.SAMTOOLS_STATS as SAMTOOLS_STATS_SINGLE_BAM {
            input:                 
                in_bam_file=MERGE_BAM.merged_bam_file            
        }
    }

    if (!OUTPUT_SINGLE_BAM) {
        Array[File] output_calling_bam_files = calling_bam
        Array[File] output_calling_bam_index_files = calling_bam_index
    }    

    if (OUTPUT_GAF){
        call utils.MERGE_GAF {
            input:
                in_sample_name=SAMPLE_NAME,
                in_gaf_chunk_files=gaf_chunks
        }        
    }    
    
################################################################
# LIFT GENOME ANNOTATIONS #
################################################################

    File position_report = select_first([LIFTOVER_POSITION_REPORT,[]])
    File liftover_vcf_script = select_first([LIFTOVER_VCF_SCRIPT,[]])
    File liftover_reference_index = select_first([LIFTOVER_LINEAR_REFERENCE_INDEX,[]])
    
    if (LIFTOVER_VCF) {
        call utils.LIFTOVER_VCF_FILE as LIFTOVER_DEEPVARIANT {
            input:
                in_vcf_file=CONCATCLIPPEDVCFCHUNKS_DEEPVARIANT.output_merged_vcf,
                in_position_report=position_report,
                script=liftover_vcf_script
        }
        call utils.REHEADER_VCF as REHEADERVCF_DEEPVARIANT {
            input:
                in_vcf_file=LIFTOVER_DEEPVARIANT.vcf_liftover,
                in_reference_index_file=liftover_reference_index
        }             
    }

################################################################
# FILTER VARIANTS #
################################################################
      
    File variantcaller_DeepVariant_VCF_output = select_first([REHEADERVCF_DEEPVARIANT.reheader_vcf, CONCATCLIPPEDVCFCHUNKS_DEEPVARIANT.output_merged_vcf])
    File variantcaller_DeepVariant_VCF_output_index = select_first([REHEADERVCF_DEEPVARIANT.reheader_vcf_index, CONCATCLIPPEDVCFCHUNKS_DEEPVARIANT.output_merged_vcf_index])

    # Filter DeepVariant_Vcf
    call utils.VCF_FILTER as FILTER_DEEPVARIANT {
        input:
            in_vcf_file=variantcaller_DeepVariant_VCF_output,
            in_vcf_file_index=variantcaller_DeepVariant_VCF_output_index            
    }
        
################################################################
# VCF STATISTIC #
################################################################
       
    call utils.BCFTOOLS_STATS as BCFTOOLS_STATS_DEEPVARIANT {
        input:
            in_vcf_file=FILTER_DEEPVARIANT.filtered_vcf,
            in_vcf_file_index=FILTER_DEEPVARIANT.filtered_vcf_index
    }

    call utils.VCFTOOLS_SUMMARY as VCFTOOLS_SUMMARY_DEEPVARIANT {
        input:
            in_vcf_file=FILTER_DEEPVARIANT.filtered_vcf,
            in_vcf_file_index=FILTER_DEEPVARIANT.filtered_vcf_index
    }

    call utils.BCFTOOLS_STATS as BCFTOOLS_STATS_DEEPVARIANT_GVCF {
        input:
            in_vcf_file=CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT.output_merged_vcf,
            in_vcf_file_index=CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT.output_merged_vcf_index
    }

    call utils.VCFTOOLS_SUMMARY as VCFTOOLS_SUMMARY_DEEPVARIANT_GVCF {
        input:
            in_vcf_file=CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT.output_merged_vcf,
            in_vcf_file_index=CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT.output_merged_vcf_index
    }    
        
################################################################
# SNPEFF: ANNOTATIONS #
################################################################

    # Run snpEff annotation on final VCF as desired
    if (SNPEFF_ANNOTATION && defined(SNPEFF_DATABASE)) {
        call utils.SNPEFF_ANNOTATE_VCF as SNPEFF_ANNOTATE_VCF_DEEPVARIANT {
            input:
                in_vcf_file=FILTER_DEEPVARIANT.filtered_vcf,
                in_vcf_file_index=FILTER_DEEPVARIANT.filtered_vcf_index,
                in_snpeff_database=SNPEFF_DATABASE,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
        }        
    }

    if (SNPEFF_ANNOTATION_LOCAL && defined(NAME_SNPEFF_DATABASE)) {
        call utils.SNPEFF_ANNOTATE_VCF_LOCAL as SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT {
            input:
                in_vcf_file=FILTER_DEEPVARIANT.filtered_vcf,
                in_vcf_file_index=FILTER_DEEPVARIANT.filtered_vcf_index,
                in_name_snpeff_database=NAME_SNPEFF_DATABASE        
        }      
    }

################################################################
# MULTIQC: SUMMARIZE BIOINFORMATICS FINDINGS FROM MULTIPLE SAMPLES INTO ONE REPORT #
################################################################
    
    File config_multiqc = select_first([MULTIQC_CONFIG,[]])
    call utils.MULTIQC {
        input:
            input_files=select_all(
            [
                FASTP_PAIRED_READS.log,
                FASTP_PAIRED_READS.json,
                FASTP_PAIRED_READS.html,
                FASTP_SINGLE_END.log,
                FASTP_SINGLE_END.json,
                FASTP_SINGLE_END.html,
                FASTQC_PAIRED_READS.output_fastqc_html_1,
                FASTQC_PAIRED_READS.output_fastqc_report_zip_1,
                FASTQC_PAIRED_READS.output_fastqc_html_2,
                FASTQC_PAIRED_READS.output_fastqc_report_zip_2,
                FASTQC_SINGLE_END.output_fastqc_html_1,
                FASTQC_SINGLE_END.output_fastqc_report_zip_1,
                BCFTOOLS_STATS_DEEPVARIANT.bcftools_stats,
                BCFTOOLS_STATS_DEEPVARIANT_GVCF.bcftools_stats,
                SNPEFF_ANNOTATE_VCF_DEEPVARIANT.annotated_stats,
                SNPEFF_ANNOTATE_VCF_DEEPVARIANT.annotated_html,
                SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_stats,
                SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_html,
                VCFTOOLS_SUMMARY_DEEPVARIANT.TsTv_count,
                VCFTOOLS_SUMMARY_DEEPVARIANT.TsTv_qual,
                VCFTOOLS_SUMMARY_DEEPVARIANT.FILTER_summary,
                VCFTOOLS_SUMMARY_DEEPVARIANT_GVCF.TsTv_count,
                VCFTOOLS_SUMMARY_DEEPVARIANT_GVCF.TsTv_qual,
                VCFTOOLS_SUMMARY_DEEPVARIANT_GVCF.FILTER_summary,
                SAMTOOLS_STATS_SINGLE_BAM.output_samtools_stats,
            ]),
            samtools_stats=select_first([SAMTOOLS_STATS_CALLING_BAM.output_samtools_stats,[]]),
            config=config_multiqc,
            prefix=SAMPLE_NAME + ".multiqc"           
    } 

################################################################
# FINAL OUTPUTS #
################################################################    
    output {

        Array[File]? output_fastp = select_all([FASTP_PAIRED_READS.output_fastp_1_file, FASTP_PAIRED_READS.output_fastp_2_file, 
                            FASTP_PAIRED_READS.html, FASTP_SINGLE_END.output_fastp_1_file, FASTP_SINGLE_END.html])
                     
        Array[File] output_vcf_DeepVariant = select_all([CONCATCLIPPEDVCFCHUNKS_DEEPVARIANT.output_merged_vcf, CONCATCLIPPEDVCFCHUNKS_DEEPVARIANT.output_merged_vcf_index,
                            REHEADERVCF_DEEPVARIANT.reheader_vcf, REHEADERVCF_DEEPVARIANT.reheader_vcf_index,
                            FILTER_DEEPVARIANT.filtered_vcf, FILTER_DEEPVARIANT.filtered_vcf_index,
                            CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT.output_merged_vcf, CONCATCLIPPEDGVCFCHUNKS_DEEPVARIANT.output_merged_vcf_index])

        Array[File] bcftools_stats_DeepVariant = select_all([BCFTOOLS_STATS_DEEPVARIANT.bcftools_stats, BCFTOOLS_STATS_DEEPVARIANT_GVCF.bcftools_stats])
        
        Array[File]? annotated_vcf_DeepVariant = select_all([SNPEFF_ANNOTATE_VCF_DEEPVARIANT.annotated_vcf, SNPEFF_ANNOTATE_VCF_DEEPVARIANT.annotated_stats,
                            SNPEFF_ANNOTATE_VCF_DEEPVARIANT.annotated_html, SNPEFF_ANNOTATE_VCF_DEEPVARIANT.annotated_gene,
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_vcf, SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_vcf_index,
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_stats, SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_csv_report,
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_clinvar_vcf, SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_clinvar_vcf_index,
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_clinvar_csv_report, 
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_pathogenic_list, SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_pathogenic_interval,
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_pathogenic_vcf, SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_pathogenic_vcf_index,
                            SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_html, SNPEFF_ANNOTATE_VCF_LOCAL_DEEPVARIANT.annotated_gene])
        
        File? output_gaf = MERGE_GAF.output_merged_gaf
        
        Array[File]? output_single_bam = select_all([MERGE_BAM.merged_bam_file, MERGE_BAM.merged_bam_file_index])
        File? samtools_stats_single_bam = SAMTOOLS_STATS_SINGLE_BAM.output_samtools_stats

        Array[File]? output_calling_bams = output_calling_bam_files
        Array[File]? output_calling_bam_indexes = output_calling_bam_index_files
        Array[File]? samtools_stats_calling_bam = SAMTOOLS_STATS_CALLING_BAM.output_samtools_stats

        Array[File]? output_liftover_bams = output_calling_bam_liftover_files
        Array[File]? output_liftover_bam_indexes = output_calling_bam_liftover_index_files
        
        File multiqc = MULTIQC.multiqc_report      
    }
}
