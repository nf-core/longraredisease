#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    longraredisease: Comprehensive Nanopore Rare Disease Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import nf-schema function
include { samplesheetToList } from 'plugin/nf-schema'


// Data preprocessing subworkflows
include { BAM_STATS_SAMTOOLS                 } from '../subworkflows/nf-core/bam_stats_samtools/main.nf'
include { SAMTOOLS_INDEX                     } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                     } from '../modules/nf-core/samtools/faidx/main.nf'
include { bam2fastq                          } from '../subworkflows/local/bam2fastq.nf'
include { align                              } from '../subworkflows/local/align.nf'
include { CAT_FASTQ                          } from '../modules/nf-core/cat/fastq/main.nf'
include { NANOPLOT as NANOPLOT_QC            } from '../modules/nf-core/nanoplot/main'
include { CREATE_PEDIGREE_FILE               } from '../modules/local/create_ped_file/main.nf'

// Coverage analysis subworkflows
include { mosdepth                           } from '../subworkflows/local/mosdepth.nf'
include { multiqc_mosdepth                   } from '../subworkflows/local/multiqc_mosdepth.nf'

// Trio analysis - rtg format reference file
include { RTG_FORMAT_REF                     } from '../modules/local/rtg/format_ref/main.nf'

// Methylation calling
include { methyl                             } from '../subworkflows/local/methyl.nf'

// SNV/indel calling
include { call_snv                           } from '../subworkflows/local/call_snv'
// annotate snv - future releases

// Haplotag BAM
include { SNIFFLES as SNIFFLES_UNPHASED      } from '../modules/nf-core/sniffles/main.nf'
include { longphase_variants                  } from '../subworkflows/local/longphase_variants.nf'
include { haplotag_bam                       } from '../subworkflows/local/haplotag_bam.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAPLOTAG } from '../modules/nf-core/samtools/index/main'

// SV calling
include { call_sv                            } from '../subworkflows/local/call_sv.nf'
include { SNIFFLES_GENERATE_PLOTS            } from '../modules/local/sniffles/generate_plots/main.nf'
include { filter_sv as filter_sv_sniffles    } from '../subworkflows/local/filter_sv'

// Merge SV - multiple callers
include { filter_sv as filter_sv_svim        } from '../subworkflows/local/filter_sv'
include { filter_sv as filter_sv_cutesv      } from '../subworkflows/local/filter_sv'
include { GUNZIP as GUNZIP_SVIM              } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_CUTESV            } from '../modules/nf-core/gunzip/main.nf'
include { merge_sv                           } from '../subworkflows/local/merge_sv.nf'

// Annotate and prioritize variants
include { annotsv_db                         } from '../subworkflows/local/annotsv_db.nf'
include { annotate_sv                        } from '../subworkflows/local/annotate_sv.nf'
include { SVANNA_PRIORITIZE                  } from '../modules/local/svanna/main.nf'

// SV calling for trios
include { sniffles_trio                      } from '../subworkflows/local/sniffles_trio.nf'
include { rtg_compare_sv                     } from '../subworkflows/local/rtg_compare_sv.nf'
include { annotate_sv as annotate_mendelian_sv  } from '../subworkflows/local/annotate_sv.nf'
include { annotate_sv as annotate_denovo_sv  } from '../subworkflows/local/annotate_sv.nf'

// SNV calling for trios
include { joint_genotype_snv                 } from '../subworkflows/local/joint_genotype_snv.nf'
include { rtg_compare_snv                    } from '../subworkflows/local/rtg_compare_snv.nf'

// STR analysis subworkflow
include { call_str                           } from '../subworkflows/local/call_str.nf'
include { annotate_str                       } from '../subworkflows/local/annotate_str.nf'

// CNV calling subworkflows
include { call_spectre_cnv                   } from '../subworkflows/local/call_spectre_cnv.nf'
include { call_hificnv                       } from '../subworkflows/local/call_hificnv.nf'

// VCF processing subworkflows
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_longraredisease_pipeline'


workflow longraredisease {

    main:

    // Convert samplesheet to list and create channel using nf-schema
    def samplesheet_data = samplesheetToList(params.input, "assets/schema_input.json")

    ch_samplesheet = Channel.fromList(samplesheet_data)
    .map { row ->
        // Handle the ArrayList structure from nf-schema
        if (row instanceof List) {
            def meta_map = row[0]
            def sample_id = meta_map.id ?: meta_map.toString()
            def meta = [id: sample_id]
            def data = [
                id: sample_id,
                file_path: row[1],
                hpo_terms: row[2] ?: null,
                sex: row[3] ?: 0,
                phenotype: row[4] ?: 0,
                family_id: row[5] ?: null,
                maternal_id: row[6] ?: "0",
                paternal_id: row[7] ?: "0"
                ]
            return [meta, data]
        } else {
            error "Unexpected row type: ${row.getClass()}"
        }
    }

    if (params.trio_analysis) {

        pedigree_input = ch_samplesheet
        .map { meta, data ->
        // Extract family_id from data, not from sample
        [data.family_id, data]
        }
        .groupTuple(size:3)  // Groups by first element (family_id) does not have a size parameter (blocking operation)
        // 3 samples size 3 - as soon as 3 samples finish then it can continue group key with group tuple (compute teh size when you don't knwo the amount of samples)
        .map { family_id, family_samples ->
        def family_meta = [id: family_id]
        [family_meta, family_samples]
        }

        CREATE_PEDIGREE_FILE(pedigree_input)

        }
/*
=======================================================================================
                                REFERENCE FILES SETUP
=======================================================================================
*/
    // Initialize versions channel
    ch_versions = Channel.empty()

    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }
        .first()

    // Generate FAI index
    SAMTOOLS_FAIDX(ch_fasta, [[:], []], true)

    ch_fai = SAMTOOLS_FAIDX.out.fai

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Create combined FASTA+FAI channel by joining

    ch_fasta_fai = ch_fasta
        .join(ch_fai, by: 0)
        .map { meta, fasta, fai -> tuple(meta, fasta, fai) }
        .first()

    if (params.trio_analysis){

        RTG_FORMAT_REF(ch_fasta)  // Run format_ref early to get SDF for trio comparison

        // Extract unique family IDs from samplesheet
        ch_family_ids = ch_samplesheet
            .map { meta, data -> data.family_id }
            .unique()

        // Replicate SDF with family-specific metadata
        ch_sdf = RTG_FORMAT_REF.out.sdf
            .map { meta, sdf -> sdf }  // Extract just the SDF path
            .combine(ch_family_ids)     // Combine with each family_id
            .map { sdf, family_id -> [[id: family_id], sdf] }


            }

    // Tandem repeat file for Sniffles (only if SV calling is enabled)
    if (params.sv) {

        ch_trf = Channel
            .fromPath(params.sniffles_tandem_file, checkIfExists: true)
            .map { bed -> tuple([id: "trf"], bed) }
            .first()

            }

            else {

                ch_trf = Channel.empty()

                }


/*
=======================================================================================
                                DATA PREPROCESSING PIPELINE
=======================================================================================
*/
    if (params.input_type == 'fastq') {
        /*
        ================================================================================
                            FASTQ ALIGNMENT WORKFLOW
        ================================================================================
        */
        // Collect FASTQ files
        ch_fastq_files = ch_samplesheet
        .map { meta, data ->
            def fastq = file(data.file_path)

            if (fastq.isFile() && (fastq.name.endsWith('.fastq.gz') || fastq.name.endsWith('.fq.gz'))) {
                // Single FASTQ file case
                return [meta, [fastq]]
            } else if (fastq.isDirectory()) {
                // Directory with multiple FASTQ files
                def fastq_files = fastq.listFiles().findAll {
                    it.name.endsWith('.fastq.gz') || it.name.endsWith('.fq.gz')
                }

                if (fastq_files.isEmpty()) {
                    error "No FASTQ files found in directory: ${data.fastq} for sample ${meta.id}"
                }

                return [meta, fastq_files]
            } else {
                error "Invalid FASTQ input for sample ${meta.id}: ${data.fastq}"
            }
        }

        ch_fastq_files.branch { meta, files ->
        single: files.size() == 1
            return [meta, files[0]]  // Extract single file from list
            multiple: files.size() > 1
            return [meta + [single_end: true], files]  // Keep as list for CAT_FASTQ
            }.set { fastq_branched }


        // Prepare input for nanoplot from FASTQ
        CAT_FASTQ(
            fastq_branched.multiple.map { meta, fastq_list ->
                [meta + [single_end: true], fastq_list]
            }
        )

        ch_processed_fastq = fastq_branched.single
        .mix(CAT_FASTQ.out.reads)

        // Align FASTQ reads to reference genome using minimap2
        align (
            ch_fasta,
            ch_processed_fastq,
            params.winnowmap_kmers,
            params.filter_targets
        )

        ch_versions = ch_versions.mix(align.out.versions)

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = align.out.bam
        ch_final_sorted_bai = align.out.bai

        ch_nanoplot = ch_processed_fastq
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    }

    else if (params.input_type == 'ubam') {
        /*
        ================================================================================
                            ALIGNMENT WORKFLOW (UNALIGNED INPUT)
        ================================================================================
        */


        // Collect unaligned BAM files
        ch_bam_files = ch_samplesheet
            .map { meta, data ->
                def bam_input = data.file_path

                if (!bam_input) {
            error "No BAM input provided for sample ${meta.id}"
                }

                def bam_path = file(bam_input)

                if (bam_path.isFile() && bam_path.name.endsWith('.bam')) {
                    // Single BAM file case
                    return [meta + [is_multiple: false], bam_path]
                } else if (bam_path.isDirectory()) {
                    // Directory with multiple BAM files case
                    def bam_files = bam_path.listFiles().findAll { it.name.endsWith('.bam') }

                    if (bam_files.isEmpty()) {
                        error "No BAM files found for sample ${meta.id} in directory: ${bam_input}"
                    }

                    return [meta + [is_multiple: bam_files.size() > 1], bam_files]
                } else {
                    error "Invalid BAM input for sample ${meta.id}: ${bam_input} (not a file or directory)"
                }
            }

        // Convert BAM to FASTQ
        bam2fastq (
            ch_bam_files,
            [[:], []],
            [[:], []],
            [[:], []]
        )

        ch_versions = ch_versions.mix(bam2fastq.out.versions)
        // Align FASTQ reads to reference genome using minimap2
        align (
            ch_fasta,
            bam2fastq.out.other,
            params.winnowmap_kmers,
            params.filter_targets
        )

        ch_versions = ch_versions.mix(align.out.versions)

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = align.out.bam
        .map { meta, bam ->
        def clean_meta = [id: meta.id]
        [clean_meta, bam]
        }

        ch_final_sorted_bai = align.out.bai
        .map { meta, bai ->
        def clean_meta = [id: meta.id]
        [clean_meta, bai]
        }


        // Prepare input for nanoplot from FASTQ
        ch_nanoplot = bam2fastq.out.other
            .map { meta, fastq_file ->
                tuple(meta, fastq_file)
            }

    } else if (params.input_type == 'bam') {
        /*
        ================================================================================
                            ALIGNED INPUT WORKFLOW (ALIGNED BAM INPUT)
        ================================================================================
        */

        // For aligned BAM input
        ch_aligned_input = ch_samplesheet
            .map { meta, data ->
                def bam_file = file(data.file_path, checkIfExists: true)
                def bai_file = file("${data.file_path}.bai", checkIfExists: true)
                return [meta, bam_file, bai_file]
            }

        // Use this single channel for all downstream processes
        ch_final_sorted_bam = ch_aligned_input.map { meta, bam, bai -> [meta, bam] }
        ch_final_sorted_bai = ch_aligned_input.map { meta, bam, bai -> [meta, bai] }

        // For nanoplot, we'll skip it (no FASTQ available)
        ch_nanoplot = ch_final_sorted_bam
    }

        // Prepare input channel with BAM, BAI, and optional BED file for coverage analysis
    ch_input_bam_bai_bed = ch_final_sorted_bam
        .join(ch_final_sorted_bai, by: 0)
        .map { meta, bam, bai ->
            def bed = params.targets_bed ? file(params.targets_bed) : []
            tuple(meta, bam, bai, bed)
        }

    // Prepare simplified BAM input channel for variant calling and methylation calling
    ch_input_bam = ch_final_sorted_bam
        .join(ch_final_sorted_bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }


/*
=======================================================================================
                                METHYLATION ANALYSIS
=======================================================================================
*/

    if (params.methyl) {
        // Use workflow-generated BAM for methylation analysis
        ch_methyl_input = ch_input_bam

        methyl(
            ch_methyl_input,
            ch_fasta_fai,
            [[:], []]
            )

            ch_versions = ch_versions.mix(methyl.out.versions)

            }

/*
=======================================================================================
                                QC
=======================================================================================
*/


    if (params.generate_bam_stats) {
        // Run BAM statistics using samtools
        BAM_STATS_SAMTOOLS(
            ch_input_bam,
            ch_fasta
        )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    }

    // Run nanoplot (only if we have FASTQ data from alignment workflow)
    if (params.qc) {

        NANOPLOT_QC(
            ch_nanoplot
        )
        ch_versions = ch_versions.mix(NANOPLOT_QC.out.versions)
    }


/*
=======================================================================================
                                Coverage analysis
=======================================================================================
*/
    // Run mosdepth when needed
    if (params.downsample_sv || params.generate_coverage_report || params.cnv_spectre || params.cnv_hificnv) {

        mosdepth(
            ch_input_bam_bai_bed,
            [[:], []]
        )
        ch_versions = ch_versions.mix(mosdepth.out.versions)

        // Combine all mosdepth outputs per sample, preserving metadata

        ch_mosdepth = mosdepth.out.global_txt
        .join(mosdepth.out.summary_txt)
        .join(mosdepth.out.regions_txt)
        .map { meta, file1, file2, file3 ->
        [meta, [file1, file2, file3]]  // Combine files into a single list
        }

        multiqc_mosdepth (
            ch_mosdepth  // Pass [meta, [files]] tuples
        )

        ch_versions = ch_versions.mix(multiqc_mosdepth.out.versions)

        }


/*
=======================================================================================
                                CALL SNV/INDEL
=======================================================================================
*/

    if (params.snv || params.haplotag_bam) {

        call_snv (
            ch_input_bam,
            ch_fasta,
            ch_fai,
            params.run_deepvariant,
            ch_input_bam_bai_bed
        )

        ch_versions = ch_versions.mix(call_snv.out.versions)

        ch_snv_vcf = call_snv.out.vcf
        ch_snv_phased_vcf = call_snv.out.phased_vcf


        }

/*
======================================================================================================
                            HAPLOTAG BAM - required for SV, STR and HIFICNV so automatically turned on
======================================================================================================
*/

    if (params.haplotag_bam){

        SNIFFLES_UNPHASED(
            ch_input_bam,
            ch_fasta,
            ch_trf,
            params.vcf_output,
            params.snf_output
        )

        longphase_variants(
            ch_input_bam,
            ch_snv_vcf,
            SNIFFLES_UNPHASED.out.vcf,
            ch_fasta,
            ch_fai
        )

        haplotag_bam(
        ch_input_bam,
        longphase_variants.out.snv_vcf,
        longphase_variants.out.sv_vcf,
        ch_fasta,
        ch_fai
        )

        SAMTOOLS_INDEX_HAPLOTAG(haplotag_bam.out.bam)

        ch_input_bam = haplotag_bam.out.bam
        .join(SAMTOOLS_INDEX_HAPLOTAG.out.bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }

        ch_versions = ch_versions.mix(SNIFFLES_UNPHASED.out.versions)
        ch_versions = ch_versions.mix(longphase_variants.out.versions)
        ch_versions = ch_versions.mix(haplotag_bam.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HAPLOTAG.out.versions)

        }


/*
=======================================================================================
                                CALL SV
=======================================================================================
*/
    if (params.sv){

        call_sv(
            ch_input_bam,
            ch_fasta,
            ch_trf,
            params.vcf_output,
            params.snf_output,
            params.merge_sv,
        )
        ch_versions = ch_versions.mix(call_sv.out.versions)

        // Initialize SV VCF channels WITH CALLER INFO in metadata

        ch_final_sv_vcf = call_sv.out.sniffles_vcf
        ch_svim_vcf = call_sv.out.svim_vcf
        ch_cutesv_vcf = call_sv.out.cutesv_vcf



    }


    /*
    ================================================================================
                            OPTIONAL SV FILTERING
    ================================================================================
    */

    if (params.filter_pass_sv) {
        // Only filter the callers that were actually run
        if (params.sv) {
            filter_sv_sniffles(
                call_sv.out.sniffles_vcf_tbi
                    .filter { meta, vcf, tbi -> vcf != null }
                    .map { meta, vcf, tbi -> [meta + [caller: 'sniffles'], vcf, tbi] },
                params.coverage_bed,
                params.downsample_sv,
                mosdepth.out.summary_txt,
                mosdepth.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support,
                params.min_read_support_limit
            )
            ch_final_sv_vcf = filter_sv_sniffles.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
            ch_versions = ch_versions.mix(filter_sv_sniffles.out.versions)
        }

        if (params.sv && params.merge_sv) {
            filter_sv_svim(
                call_sv.out.svim_vcf_tbi
                    .filter { meta, vcf, tbi -> vcf != null }
                    .map { meta, vcf, tbi -> [meta + [caller: 'svim'], vcf, tbi] },
                params.coverage_bed,
                params.downsample_sv,
                mosdepth.out.summary_txt,
                mosdepth.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support,
                params.min_read_support_limit
            )
            ch_svim_vcf = filter_sv_svim.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
            ch_versions = ch_versions.mix(filter_sv_svim.out.versions)
        }

        if (params.sv && params.merge_sv) {
            filter_sv_cutesv(
                call_sv.out.cutesv_vcf_tbi
                    .filter { meta, vcf, tbi -> vcf != null }
                    .map { meta, vcf, tbi -> [meta + [caller: 'cutesv'], vcf, tbi] },
                params.coverage_bed,
                params.downsample_sv,
                mosdepth.out.summary_txt,
                mosdepth.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support,
                params.min_read_support_limit
            )
            ch_cutesv_vcf = filter_sv_cutesv.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
            ch_versions = ch_versions.mix(filter_sv_cutesv.out.versions)
        }
    }

    /*
    ================================================================================
                            MERGE SV - optional
    ================================================================================
    */

    // Gunzip VCFs for Jasmine (requires uncompressed input)
    if (params.merge_sv){
        GUNZIP_SVIM(ch_svim_vcf)
        GUNZIP_CUTESV(ch_cutesv_vcf)
        ch_versions = ch_versions.mix(GUNZIP_SVIM.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_CUTESV.out.versions)

        // Prepare input for JASMINESV - group all uncompressed VCFs by sample
        jasmine_input_ch = call_sv.out.sniffles_unzipped_vcf
            .map { meta, vcf -> [[id: meta.id], vcf] }
            .join(
                GUNZIP_SVIM.out.gunzip.map { meta, vcf -> [[id: meta.id], vcf] },
                by: 0
            )
            .join(
                GUNZIP_CUTESV.out.gunzip.map { meta, vcf -> [[id: meta.id], vcf] },
                by: 0
            )
            .map { sample_key, sniffles_vcf, svim_vcf, cutesv_vcf ->
                [sample_key, [sniffles_vcf, svim_vcf, cutesv_vcf]]
            }
            .join(
                ch_input_bam.map { meta, bam, bai -> [[id: meta.id], bam, bai] },
                by: 0
            )
            .map { sample_key, vcfs, bam, bai ->
                def clean_meta = [id: sample_key.id]
                [clean_meta, vcfs, bam, bai, []]  // [meta, vcfs, bam, bai, sample_dists]
            }

        // Run JASMINESV merging
        merge_sv(
            jasmine_input_ch,
            ch_fasta,
            ch_fai,
            []
        )
        ch_versions = ch_versions.mix(merge_sv.out.versions)

        // Set final SV VCF to merged result
        ch_final_sv_vcf = merge_sv.out.vcf
            .map { meta, vcf -> [meta + [caller: 'merged'], vcf] }
    }


    /*
    ================================================================================
                            Annotate SV vcf
    ================================================================================
    */

    if (params.sv && params.annotate_sv){

        annotsv_db(params.annotsv_annotations)

        ch_hpo_terms = ch_samplesheet.map { meta, data ->
        [meta, data.hpo_terms]
        }

        // Call the subworkflow
        annotate_sv(
            ch_final_sv_vcf,       // [[id:test, caller:sniffles], vcf, index]
            ch_hpo_terms,        // [[id:test], "HP:0001249,HP:0001250"] or [[id:test], ""]
            ch_snv_vcf,          // [[id:test], snv_vcf, snv_index]
            annotsv_db.out.db,
            [],
            [],
            []
            )

            ch_versions = ch_versions.mix(annotate_sv.out.versions)

            }

    /*
    ================================================================================
                            SV ANNOTATION WITH SVANNA
    ================================================================================
    */

    if (params.sv && params.run_svanna) {
        // Filter samplesheet to only include samples with HPO terms
        ch_samplesheet_with_hpo = ch_samplesheet
            .filter { meta, data ->
                data.hpo_terms && data.hpo_terms.trim() != ""
            }

        ch_hpo_terms = ch_samplesheet_with_hpo.map { meta, data ->
            [meta, data.hpo_terms]
        }

        // Prepare VCF for annotation with HPO terms
        ch_sv_vcf_for_annotation = ch_final_sv_vcf
            .map { meta, vcf -> [meta.id, vcf, meta.caller] }
            .join(ch_hpo_terms.map { meta, hpo -> [meta.id, hpo] }, by: 0)
            .map { sample_id, vcf, caller, hpo_terms ->
                def meta = [id: sample_id, caller: caller]
                [meta, vcf, hpo_terms]
            }

        // Set up Svanna database
        ch_svanna_db = Channel
            .fromPath(params.svanna_db, checkIfExists: true)
            .first()


        SVANNA_PRIORITIZE(
            ch_sv_vcf_for_annotation.map { meta, vcf, hpo_terms -> [meta, vcf] },
            ch_svanna_db,
            ch_sv_vcf_for_annotation.map { meta, vcf, hpo_terms -> hpo_terms }
        )

        ch_versions = ch_versions.mix(SVANNA_PRIORITIZE.out.versions)
    }

/*
=======================================================================================
                                Trio analysis
=======================================================================================
*/

    if (params.sv && params.trio_analysis) {
        sniffles_trio(call_sv.out.sniffles_snf,
        ch_samplesheet,
        ch_fasta)


    ch_trio_sv_vcf = sniffles_trio.out.vcf
        .map { meta, vcf -> [meta + [variant_type: 'sv'], vcf] }

    rtg_compare_sv(
            ch_sdf,
            ch_trio_sv_vcf,
            CREATE_PEDIGREE_FILE.out.ped
                .map { meta, ped -> [meta, ped] },
                params.run_mendelian,
                params.run_denovo
            )


            }

    if (params.snv && params.trio_analysis) {

        ch_gvcf = call_snv.out.gvcf

        joint_genotype_snv(
            ch_gvcf,
            ch_samplesheet,
            [[:], []]  // ch_bed (empty)
            )

        ch_versions = ch_versions.mix(joint_genotype_snv.out.versions)

        ch_trio_snv_vcf = joint_genotype_snv.out.vcf
        .map { meta, vcf -> [meta + [variant_type: 'snv'], vcf] }

        rtg_compare_snv(
            ch_sdf,
            ch_trio_snv_vcf,
            CREATE_PEDIGREE_FILE.out.ped
            .map { meta, ped -> [meta, ped] },
            params.run_mendelian,
            params.run_denovo
        )


        }

    // annotate trios - future release


/*
===============================================================================
                        SHORT TANDEM REPEAT ANALYSIS
================================================================================
*/

    if (params.str) {
        call_str (
            ch_input_bam,
            ch_fasta,
            params.str_bed_file
        )

        ch_variant_catalogue = channel.fromPath(params.variant_catalogue)
        .map { file -> [ [id: 'variant_catalog'], file ] }
        .first()

        annotate_str(
            call_str.out.vcf,
            ch_variant_catalogue
        )

        ch_versions = ch_versions.mix(call_str.out.versions)

        }

/*
=======================================================================================
                        COPY NUMBER VARIANT CALLING
=======================================================================================
*/

    if (params.sequencing_platform== 'pacbio' && params.cnv || params.sequencing_platform== 'hifi' && params.cnv || params.filter_targets && params.cnv) {

        ch_bam_bai_maf = ch_input_bam
        .join(ch_snv_phased_vcf.map { meta, vcf -> [[id: meta.id], vcf] }, by: 0)
        .map { meta, bam, bai, maf -> [meta, bam, bai, maf] }

        // Create channels for exclude bed (optional)
        ch_exclude = params.hificnv_exclude_bed
         ? channel.of([[id: 'exclude'], file(params.hificnv_exclude_bed, checkIfExists: true)]).first()
         : channel.of([[id: 'exclude'], []]).first()


        // Create channels for expected CN bed (optional)
        ch_expected_cn = params.hificnv_expected_cn_bed
        ? channel.of([[id: 'expected_cn'], file(params.hificnv_expected_cn_bed, checkIfExists: true)]).first()
        : channel.of([[id: 'expected_cn'], []]).first()


        call_hificnv(
            ch_bam_bai_maf,
            ch_fasta,
            ch_exclude,
            ch_expected_cn
        )

        ch_cnv_vcf = call_hificnv.out.vcf
        ch_versions = ch_versions.mix(call_hificnv.out.versions)
    }

    if (params.sequencing_platform == 'ont' && params.cnv && !params.filter_targets) {

        if (params.use_test_data) {

            // Use test data as does not accept filtered vcfs - needs whole genome to work

            ch_test_meta = channel.of([id: 'test'])

            ch_test_summary = ch_test_meta.map { meta -> [meta, file(params.spectre_test_summary_txt)]}
            ch_test_regions_bed = ch_test_meta.map { meta -> [meta, file(params.spectre_test_regions_bed)]}
            ch_test_regions_csi = ch_test_meta.map { meta -> [meta, file(params.spectre_test_regions_csi)]}
            ch_test_vcf = ch_test_meta.map { meta -> [meta, file(params.spectre_test_clair3_vcf)]}
            ch_test_fasta = ch_test_meta.map { meta -> [meta, file(params.spectre_test_fasta_file)]}

            call_spectre_cnv(
                ch_test_summary,
                ch_test_regions_bed,
                ch_test_regions_csi,
                ch_test_vcf,
                ch_test_fasta,
                params.spectre_metadata,
                params.spectre_blacklist,
                1000
            )

            ch_cnv_vcf = call_spectre_cnv.out.vcf
            ch_versions = ch_versions.mix(call_spectre_cnv.out.versions)
            }

        else {

            call_spectre_cnv(
                mosdepth.out.summary_txt,
                mosdepth.out.regions_bed,
                mosdepth.out.regions_csi,
                ch_snv_vcf,
                ch_fasta,
                params.spectre_metadata,
                params.spectre_blacklist,
                params.spectre_bin_size ?: 1000
                )

                ch_cnv_vcf = call_spectre_cnv.out.vcf
                ch_versions = ch_versions.mix(call_spectre_cnv.out.versions)

                }
                }

// Collect all versions and generate YAML
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_longraredisease_software_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_collated_versions



}
