version 1.0

import "../tasks/bioinfo_gatk_utils.wdl" as utils
import "../tasks/gam_gaf_utils.wdl" as gautils
import "../tasks/vg_map_hts.wdl" as map
import "../tasks/deepvariant.wdl" as dv

workflow GiraffeDeepVariantGATK {

    meta {
        description: "## Giraffe-DeepVariant-GATK workflow \n The full workflow to go from sequencing reads (FASTQs, CRAM) to small variant calls (VCF). Reads are mapped to a pangenome with vg giraffe and pre-processed (e.g. indel realignment). DeepVariant then calls small variants. More information at [https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-deepvariant-workflow](https://github.com/vgteam/vg_wdl/tree/gbz#giraffe-deepvariant-workflow)."
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
        OTHER_MAKEEXAMPLES_ARG: "Additional arguments for the make_examples step of DeepVariant"
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
        Boolean OUTPUT_SINGLE_BAM = false
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
        Int SPLIT_READ_CORES = 8
        Int MAP_CORES = 16
        Int MAP_MEM = 120
        Int CALL_CORES = 8
        Int CALL_MEM = 50
    }

    if (defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
	    call utils.convertCRAMtoFASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_paired_reads=PAIRED_READS,
            in_cores=SPLIT_READ_CORES
	    }
    }
    
    if (TRIM_FASTQ && PAIRED_READS) {
        call utils.runFASTP_PAIRED_READS {
            input:
                in_sample_name = SAMPLE_NAME,
                in_read_1_fastq = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file]),
                in_read_2_fastq = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file]),
                in_call_mem = CALL_MEM
        }
    } 
    
    if (TRIM_FASTQ && !PAIRED_READS) {
        call utils.runFASTP_SINGLE_END {
            input:
                in_sample_name = SAMPLE_NAME,
                in_read_1_fastq = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file]),
                in_call_mem = CALL_MEM
        }
    }     

    File read_1_file = select_first([runFASTP_PAIRED_READS.output_fastp_1_file, runFASTP_SINGLE_END.output_fastp_1_file, INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file])
    
    # Split input reads into chunks for parallelized mapping
    call utils.splitReads as firstReadPair {
        input:
            in_read_file=read_1_file,
            in_pair_id="1",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
    }
    
    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from GBZ file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call map.extractSubsetPathNames {
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
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractSubsetPathNames.output_path_list_file, written_path_names_file])
    
    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph, we generate it ourselves, from the graph.
    if (!defined(REFERENCE_FILE)) {
        call map.extractReference {
            input:
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_prefix_to_strip=REFERENCE_PREFIX,
            in_extract_mem=MAP_MEM
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.indexReference {
            input:
                in_reference_file=reference_file
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])
    
    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################

    if (PAIRED_READS) {
        File read_2_file = select_first([runFASTP_PAIRED_READS.output_fastp_2_file, INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file])
        call utils.splitReads as secondReadPair {
            input:
            in_read_file=read_2_file,
            in_pair_id="2",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES
        }
        Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
        scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
            call map.runVGGIRAFFE as runVGGIRAFFEpe {
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
            call map.runVGGIRAFFE as runVGGIRAFFEse {
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

    Array[File] gaf_chunks = select_first([runVGGIRAFFEpe.chunk_gaf_file, runVGGIRAFFEse.chunk_gaf_file])
    scatter (gaf_file in gaf_chunks) {
        call gautils.surjectGAFtoBAM {
            input:
            in_gaf_file=gaf_file,
            in_gbz_file=GBZ_FILE,
            in_path_list_file=pipeline_path_list_file,
            in_sample_name=SAMPLE_NAME,
            in_max_fragment_length=MAX_FRAGMENT_LENGTH,
            in_paired_reads=PAIRED_READS,
            mem_gb=MAP_MEM
        }

        call utils.sortBAM {
            input:
            in_bam_file=surjectGAFtoBAM.output_bam_file,
            in_ref_dict=reference_dict_file,
            in_prefix_to_strip=REFERENCE_PREFIX
        }
    }

    call utils.mergeAlignmentBAMChunks {
        input:
        in_sample_name=SAMPLE_NAME,
        in_alignment_bam_chunk_files=sortBAM.sorted_bam
    }

    # Split merged alignment by contigs list
    call utils.splitBAMbyPath {
        input:
        in_sample_name=SAMPLE_NAME,
        in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
        in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
        in_path_list_file=pipeline_path_list_file,
        in_prefix_to_strip=REFERENCE_PREFIX
    }

    ################################################################
    # Call variants in each contig #
    ################################################################

    scatter (bam_and_index_for_path in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
        ## Evantually shift and realign reads
        if (LEFTALIGN_BAM) {
            # Just left-shift each read individually
            call utils.leftShiftBAMFile {
                input:
                in_bam_file=bam_and_index_for_path.left,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file
            }
        }
        if (REALIGN_INDELS) {
            File forrealign_bam = select_first([leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
            File forrealign_index = select_first([leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])
            # Do indel realignment
            call utils.prepareRealignTargets {
                input:
                in_bam_file=forrealign_bam,
                in_bam_index_file=forrealign_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_reference_dict_file=reference_dict_file,
                in_expansion_bases=REALIGNMENT_EXPANSION_BASES
            }
            call utils.runAbraRealigner {
                input:
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_target_bed_file=prepareRealignTargets.output_target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    # If the user has set a very low memory for mapping, don't use more for realignment
                    mem_gb=if MAP_MEM < 40 then MAP_MEM else 40
            }
        }
        File calling_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.output_bam_file, bam_and_index_for_path.left])
        File calling_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, leftShiftBAMFile.output_bam_index_file, bam_and_index_for_path.right])

        ## DeepVariant calling
        call dv.runDeepVariantMakeExamples {
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
        call dv.runDeepVariantCallVariants {
            input:
                in_sample_name=SAMPLE_NAME,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_examples_file=runDeepVariantMakeExamples.examples_file,
                in_nonvariant_site_tf_file=runDeepVariantMakeExamples.nonvariant_site_tf_file,
                in_model_meta_file=DV_MODEL_META,
                in_model_index_file=DV_MODEL_INDEX,
                in_model_data_file=DV_MODEL_DATA,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
        }

        ## GATKHaplotypeCaller calling
        call utils.runGATKHaplotypeCaller {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=calling_bam,
                in_bam_file_index=calling_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_reference_dict_file=reference_dict_file,
                in_call_cores=CALL_CORES,
                in_call_mem=CALL_MEM
            }
    }

    # Merge distributed variant called VCFs
    call utils.concatClippedVCFChunks as concatClippedVCFChunks_DeepVariant {
        input:
            in_sample_name=SAMPLE_NAME,
            in_tools="deepvariant",
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_vcf_file
    }
    call utils.concatClippedVCFChunks as concatClippedGVCFChunks_DeepVariant {
        input:
            in_sample_name=SAMPLE_NAME,
            in_tools="deepvariant",
            in_clipped_vcf_chunk_files=runDeepVariantCallVariants.output_gvcf_file
    }    
    call utils.concatClippedVCFChunks as concatClippedVCFChunks_GATK {
        input:
            in_sample_name=SAMPLE_NAME,
            in_tools="GATK",
            in_clipped_vcf_chunk_files=runGATKHaplotypeCaller.genotyped_vcf               
    }

    ################################################################
    # Filer variants #
    ################################################################

    # Filter DeepVariant_Vcf
    call utils.runBCFfilter as runDeepVariantfiltering {
        input:
            in_sample_name=SAMPLE_NAME,
            in_vcf_file=concatClippedVCFChunks_DeepVariant.output_merged_vcf,
            in_vcf_file_index=concatClippedVCFChunks_DeepVariant.output_merged_vcf_index
    }
    
    # Filter GATK_Vcf
    call utils.runGATKHardfiltering {
        input:
            in_sample_name=SAMPLE_NAME,
            in_vcf_file=concatClippedVCFChunks_GATK.output_merged_vcf,
            in_vcf_file_index=concatClippedVCFChunks_GATK.output_merged_vcf_index,
            in_reference_file=reference_file,
            in_reference_index_file=reference_index_file,
            in_reference_dict_file=reference_dict_file,
            in_call_cores=CALL_CORES,
            in_call_mem=CALL_MEM
    }
    
    if (OUTPUT_GAF){
        call gautils.mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=gaf_chunks
        }
    }

    if (OUTPUT_SINGLE_BAM){
        call utils.mergeAlignmentBAMChunks as mergeBAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=select_all(flatten([calling_bam, [splitBAMbyPath.bam_unmapped_file]]))
        }
    }

    if (!OUTPUT_SINGLE_BAM){
        Array[File] output_calling_bam_files = calling_bam
        Array[File] output_calling_bam_index_files = calling_bam_index
    }

    # FINAL OUTPUTS
    output {
        File? output_fastp1 = runFASTP_PAIRED_READS.output_fastp_1_file
        File? output_fastp2 = runFASTP_PAIRED_READS.output_fastp_2_file
        File? output_fastp_single_end = runFASTP_SINGLE_END.output_fastp_1_file
        File output_vcf_deepvariant = runDeepVariantfiltering.genotyped_filtered_vcf
        File output_vcf_index_deepvariant = runDeepVariantfiltering.genotyped_filtered_vcf_index
        File output_gvcf_deepvariant = concatClippedGVCFChunks_DeepVariant.output_merged_vcf
        File output_gvcf_index_deepvariant = concatClippedGVCFChunks_DeepVariant.output_merged_vcf_index
        File output_vcf_GATK = runGATKHardfiltering.genotyped_filtered_vcf
        File output_vcf_index_GATK = runGATKHardfiltering.genotyped_filtered_vcf_index
        File? output_gaf = mergeGAF.output_merged_gaf
        File? output_bam = mergeBAM.merged_bam_file
        File? output_bam_index = mergeBAM.merged_bam_file_index
        Array[File]? output_calling_bams = output_calling_bam_files
        Array[File]? output_calling_bam_indexes = output_calling_bam_index_files      
    }

}
