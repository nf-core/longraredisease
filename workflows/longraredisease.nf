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
include { bam2fastq_subworkflow              } from '../subworkflows/local/bam2fastq.nf'
include { alignment_subworkflow              } from '../subworkflows/local/align.nf'
include { CAT_FASTQ                          } from '../modules/nf-core/cat/fastq/main.nf'
include { NANOPLOT as NANOPLOT_QC            } from '../modules/nf-core/nanoplot/main'

// Methylation calling
include { methyl                             } from '../subworkflows/local/methyl.nf'

// Coverage analysis subworkflows
include { mosdepth                           } from '../subworkflows/local/mosdepth.nf'

// Structural variant calling subworkflows
include { call_sv                            } from '../subworkflows/local/call_sv.nf'
include { filter_sv as filter_sv_sniffles    } from '../subworkflows/local/filter_sv'
include { filter_sv as filter_sv_svim        } from '../subworkflows/local/filter_sv'
include { filter_sv as filter_sv_cutesv      } from '../subworkflows/local/filter_sv'

// SV merging and intersection filtering subworkflows
include { GUNZIP as GUNZIP_SNIFFLES          } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_SVIM              } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_CUTESV            } from '../modules/nf-core/gunzip/main.nf'
include { SV_PLOT as SV_PLOT_SNIFFLES        } from '../modules/local/generate_sv_plots/main.nf'
include { SV_PLOT as SV_PLOT_SVIM            } from '../modules/local/generate_sv_plots/main.nf'
include { SV_PLOT as SV_PLOT_CUTESV          } from '../modules/local/generate_sv_plots/main.nf'
include { merge_sv                           } from '../subworkflows/local/merge_sv.nf'
include { SVANNA_PRIORITIZE                  } from '../modules/local/SvAnna/main.nf'

// SNV calling and processing subworkflows
include { call_snv                           } from '../subworkflows/local/call_snv'
include { merge_snv_subworkflow              } from '../subworkflows/local/merge_snv.nf'

// Phasing subworkflow
include { longphase                          } from '../subworkflows/local/longphase.nf'

// CNV calling subworkflows
include { call_cnv_spectre                   } from '../subworkflows/local/call_cnv_spectre.nf'
include { call_hificnv                        } from '../subworkflows/local/call_hificnv.nf'

// STR analysis subworkflow
include { call_str                           } from '../subworkflows/local/call_str.nf'
include { annotate_str                       } from '../subworkflows/local/annotate_str.nf'

// VCF processing subworkflows
include { unify_vcf_subworkflow              } from '../subworkflows/local/unify_vcf.nf'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_longraredisease_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
                    file_path: row[1],
                    hpo_terms: row[2] ?: null,
                    sex: row[3] ?: null,
                    family_id: row[4] ?: null,
                    maternal_id: row[5] ?: null,
                    paternal_id: row[6] ?: null
                ]
                return [meta, data]
            } else {
                error "Unexpected row type: ${row.getClass()}"
            }
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

    // Tandem repeat file for Sniffles (only if SV calling is enabled)
    if (params.sv) {
        ch_trf = Channel
            .fromPath(params.sniffles_tandem_file, checkIfExists: true)
            .map { bed -> tuple([id: "trf"], bed) }
            .first()
    } else {
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
        alignment_subworkflow(
            ch_fasta,
            ch_processed_fastq,
            params.winnowmap_kmers
        )

        ch_versions = ch_versions.mix(alignment_subworkflow.out.versions)

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = alignment_subworkflow.out.bam
        ch_final_sorted_bai = alignment_subworkflow.out.bai

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
        bam2fastq_subworkflow(
            ch_bam_files,
            [[:], []],
            [[:], []],
            [[:], []]
        )

        ch_versions = ch_versions.mix(bam2fastq_subworkflow.out.versions)
        // Align FASTQ reads to reference genome using minimap2
        alignment_subworkflow(
            ch_fasta,
            bam2fastq_subworkflow.out.other,
            params.winnowmap_kmers
        )

        ch_versions = ch_versions.mix(alignment_subworkflow.out.versions)

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = alignment_subworkflow.out.bam
        .map { meta, bam ->
        def clean_meta = [id: meta.id]
        [clean_meta, bam]
        }

        ch_final_sorted_bai = alignment_subworkflow.out.bai
        .map { meta, bai ->
        def clean_meta = [id: meta.id]
        [clean_meta, bai]
        }


        // Prepare input for nanoplot from FASTQ
        ch_nanoplot = bam2fastq_subworkflow.out.other
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

/*
=======================================================================================
                                COVERAGE ANALYSIS
=======================================================================================
*/

    // Prepare input channel with BAM, BAI, and optional BED file for coverage analysis
    ch_input_bam_bai_bed = ch_final_sorted_bam
        .join(ch_final_sorted_bai, by: 0)
        .map { meta, bam, bai ->
            def bed = params.target_bed ? file(params.target_bed) : []
            tuple(meta, bam, bai, bed)
        }

    // Prepare simplified BAM input channel for variant calling and methylation calling
    ch_input_bam = ch_final_sorted_bam
        .join(ch_final_sorted_bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }

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

    // Run mosdepth when needed
    if (params.downsample_sv || params.generate_coverage) {
        mosdepth(
            ch_input_bam_bai_bed,
            [[:], []]
        )
        ch_versions = ch_versions.mix(mosdepth.out.versions)
    }

    if (params.methyl) {
            // Use workflow-generated BAM for methylation analysis
            ch_methyl_input = ch_input_bam
            ch_methyl_input.view()
        methyl(
            ch_methyl_input,
            ch_fasta_fai,
            [[:], []]
        )
        ch_versions = ch_versions.mix(methyl.out.versions)
    }

/*
=======================================================================================
                        STRUCTURAL VARIANT CALLING WORKFLOW
=======================================================================================
*/

if (params.sv) {
    /*
    ================================================================================
                            PARALLEL SV CALLER EXECUTION
    ================================================================================
    */

    // Run SV calling subworkflow
    call_sv(
        ch_input_bam,
        ch_fasta,
        ch_trf,
        params.vcf_output,
        params.snf_output
    )
    ch_versions = ch_versions.mix(call_sv.out.versions)

    // Initialize SV VCF channels WITH CALLER INFO in metadata
    ch_sniffles_vcf = call_sv.out.sniffles_vcf
    ch_svim_vcf = call_sv.out.svim_vcf
    ch_cutesv_vcf = call_sv.out.cutesv_vcf


    /*
    ================================================================================
                            OPTIONAL SV FILTERING
    ================================================================================
    */

        if (params.filter_sv_pass) {
        // Filter SVs for each caller separately - caller info automatically preserved
        filter_sv_sniffles(
        call_sv.out.sniffles_vcf_tbi.map { meta, vcf, tbi ->
            [meta + [caller: 'sniffles'], vcf, tbi]
        },
        params.target_bed,
        params.downsample_sv,
        mosdepth.out.summary_txt,
        mosdepth.out.quantized_bed,
        params.chromosome_codes,
        params.min_read_support,
        params.min_read_support_limit
    )

    filter_sv_svim(
        call_sv.out.svim_vcf_tbi.map { meta, vcf, tbi ->
            [meta + [caller: 'svim'], vcf, tbi]
        },
        params.target_bed,
        params.downsample_sv,
        mosdepth.out.summary_txt,
        mosdepth.out.quantized_bed,
        params.chromosome_codes,
        params.min_read_support,
        params.min_read_support_limit
    )

    filter_sv_cutesv(
        call_sv.out.cutesv_vcf_tbi.map { meta, vcf, tbi ->
            [meta + [caller: 'cutesv'], vcf, tbi]
        },
        params.target_bed,
        params.downsample_sv,
        mosdepth.out.summary_txt,
        mosdepth.out.quantized_bed,
        params.chromosome_codes,
        params.min_read_support,
        params.min_read_support_limit
    )

        // Update channels to use filtered results (caller info preserved automatically)
        ch_sniffles_vcf = filter_sv_sniffles.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
        ch_svim_vcf = filter_sv_svim.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
        ch_cutesv_vcf = filter_sv_cutesv.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }

        ch_versions = ch_versions.mix(filter_sv_sniffles.out.versions)
        ch_versions = ch_versions.mix(filter_sv_svim.out.versions)
        ch_versions = ch_versions.mix(filter_sv_cutesv.out.versions)
    }

    /*
    ================================================================================
                        CONDITIONAL SV MERGING OR INDIVIDUAL PROCESSING
    ================================================================================
    */

    if (params.merge_sv) {
        /*
        ========================================================================
                            SV MERGING WITH JASMINESV
        ========================================================================
        */

        // Gunzip VCFs for Jasmine (requires uncompressed input)
        GUNZIP_SNIFFLES(ch_sniffles_vcf)
        GUNZIP_SVIM(ch_svim_vcf)
        GUNZIP_CUTESV(ch_cutesv_vcf)
        ch_versions = ch_versions.mix(GUNZIP_SNIFFLES.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_SVIM.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_CUTESV.out.versions)

        // Prepare input for JASMINESV - group all uncompressed VCFs by sample
        jasmine_input_ch = GUNZIP_SNIFFLES.out.gunzip
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
        ch_sv_vcf_final = merge_sv.out.vcf
            .map { meta, vcf -> [meta + [caller: 'merged'], vcf] }

    } else {
        /*
        ========================================================================
                        INDIVIDUAL CALLER SELECTION (NO MERGING)
        ========================================================================
        */

        // Select VCF based on priority or parameter
        if (params.sv_caller == 'sniffles') {
            ch_sv_vcf_final = ch_sniffles_vcf
                .map { meta, vcf -> [meta + [caller: 'sniffles'], vcf] }
        } else if (params.sv_caller == 'svim') {
            ch_sv_vcf_final = ch_svim_vcf
                .map { meta, vcf -> [meta + [caller: 'svim'], vcf] }
        } else if (params.sv_caller == 'cutesv') {
            ch_sv_vcf_final = ch_cutesv_vcf
                .map { meta, vcf -> [meta + [caller: 'cutesv'], vcf] }
        } else {
            // Default to Sniffles if parameter not recognized
            ch_sv_vcf_final = ch_sniffles_vcf
                .map { meta, vcf -> [meta + [caller: 'sniffles'], vcf] }
        }
    }
    if (params.generate_sv_plots) {
        /*
        ================================================================================
                            SV PLOTTING
        ================================================================================
        */

        SV_PLOT_SNIFFLES(
            GUNZIP_SNIFFLES.out.gunzip.map { meta, vcf -> [meta, vcf] }
        )
        SV_PLOT_SVIM(
            GUNZIP_SVIM.out.gunzip.map { meta, vcf -> [meta, vcf] }
        )
        SV_PLOT_CUTESV(
            GUNZIP_CUTESV.out.gunzip.map { meta, vcf -> [meta, vcf] }
        )

        ch_versions = ch_versions.mix(SV_PLOT_SNIFFLES.out.versions)
    }

    /*
    ================================================================================
                            SV ANNOTATION WITH SVANNA
    ================================================================================
    */

    if (params.annotate_sv) {
        // Filter samplesheet to only include samples with HPO terms
        ch_samplesheet_with_hpo = ch_samplesheet
            .filter { meta, data ->
                data.hpo_terms && data.hpo_terms.trim() != ""
            }

        ch_hpo_terms = ch_samplesheet_with_hpo.map { meta, data ->
            [meta, data.hpo_terms]
        }

        // Prepare VCF for annotation with HPO terms
        ch_sv_vcf_for_annotation = ch_sv_vcf_final
            .map { meta, vcf -> [meta.id, vcf, meta.caller] }
            .join(ch_hpo_terms.map { meta, hpo -> [meta.id, hpo] }, by: 0)
            .map { sample_id, vcf, caller, hpo_terms ->
                def meta = [id: sample_id, caller: caller]
                [meta, vcf, hpo_terms]
            }

        // Set up SvAnna database
        ch_svanna_db = Channel
            .fromPath(params.svanna_db, checkIfExists: true)
            .first()

        // Run SvAnna prioritization
        SVANNA_PRIORITIZE(
            ch_sv_vcf_for_annotation.map { meta, vcf, hpo_terms -> [meta, vcf] },
            ch_svanna_db,
            ch_sv_vcf_for_annotation.map { meta, vcf, hpo_terms -> hpo_terms }
        )
        ch_versions = ch_versions.mix(SVANNA_PRIORITIZE.out.versions)
    }

    /*
    ================================================================================
                        SET DOWNSTREAM SV VCF CHANNEL
    ================================================================================
    */

    // Set the final SV VCF channel for downstream processes
    ch_sv_vcf_downstream = ch_sv_vcf_final
    }

    else {
    /*
    ================================================================================
                        SV CALLING DISABLED - EMPTY CHANNELS
    ================================================================================
    */

    ch_sv_vcf_downstream = Channel.empty()
    ch_sv_vcf_final = Channel.empty()

    }
/*
================================================================================
                        SINGLE NUCLEOTIDE VARIANT CALLING
================================================================================
*/

    if (params.snv) {
        // Prepare input for SNV calling
        ch_input_bam_clair3 = ch_input_bam.map { meta, bam, bai ->
            tuple(
                meta,
                bam,
                bai,
                params.clair3_model,
                [],
                params.clair3_platform
            )
        }

        // Run SNV calling
        call_snv (
            ch_input_bam_clair3,
            ch_fasta,
            ch_fai,
            params.deepvariant,
            ch_input_bam_bai_bed,
            params.filter_pass_snv
        )

        ch_versions = ch_versions.mix(call_snv.out.versions)

        ch_snv_vcf = call_snv.out.clair3_vcf
        ch_snv_tbi = call_snv.out.clair3_tbi

        if (params.merge_snv && params.deepvariant) {
            combined_vcfs = ch_snv_vcf
                .join(ch_snv_tbi, by: 0)
                .join(
                    call_snv.out.deepvariant_vcf
                        .join(call_snv.out.deepvariant_tbi, by: 0),
                    by: 0
                )
                .map { meta, clair3_vcf, clair3_tbi, deepvariant_vcf, deepvariant_tbi ->
                    [
                        meta,
                        [clair3_vcf, deepvariant_vcf],
                        [clair3_tbi, deepvariant_tbi]
                    ]
                }

            // Merge SNV VCFs
            merge_snv_subworkflow(combined_vcfs)
            ch_versions = ch_versions.mix(merge_snv_subworkflow.out.versions)
        }
    } else {
        // Create empty channels when SNV calling is disabled
        ch_snv_vcf = Channel.empty()
        ch_snv_tbi = Channel.empty()
    }

/*
=======================================================================================
                                PHASING ANALYSIS
=======================================================================================
*/

    // Run phasing with LongPhase if enabled
    if (params.phase && params.snv) {
    if (params.sv && params.phase_with_sv) {
        // Clean all metadata to just sample ID for joining
        ch_longphase_input = ch_input_bam
            .map { meta, bam, bai -> [[id: meta.id], meta, bam, bai] }
            .join(
                ch_snv_vcf.map { meta, vcf -> [[id: meta.id], vcf] },
                by: 0
            )
            .join(
                ch_sv_vcf_final.map { meta, vcf -> [[id: meta.id], vcf] },
                by: 0
            )
            .map { sample_key, original_meta, bam, bai, snv_vcf, sv_vcf ->
                tuple(original_meta, bam, bai, snv_vcf, sv_vcf, [])
            }
    } else {
        // Phasing with SNVs only
        ch_longphase_input = ch_input_bam
            .join(ch_snv_vcf, by: 0)
            .map { meta, bam, bai, snv_vcf ->
                tuple(meta, bam, bai, snv_vcf, [], [])
            }
    }

    longphase(
        ch_longphase_input,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(longphase.out.versions)

    }
/*
=======================================================================================
                        COPY NUMBER VARIANT CALLING
=======================================================================================
*/

    ch_spectre_vcf = Channel.empty()
    ch_hificnv_vcf = Channel.empty()

    if (params.cnv_spectre) {
        // Spectre CNV calling - requires SNV data unless using test data
        if (params.use_test_data) {
            // Test mode with hardcoded parameters
            ch_spectre_test_reference = ch_samplesheet
            .map { meta, data -> meta.id }  // Extract sample ID from samplesheet
            .combine(Channel.fromPath(params.spectre_test_clair3_vcf, checkIfExists: true))
            .combine(Channel.fromPath(params.spectre_test_fasta_file, checkIfExists: true))
            .map { sample_id, vcf_file, fasta ->
            def meta = [id: sample_id]
            tuple(meta, fasta)
            }

            call_cnv_spectre(
                params.spectre_test_mosdepth,
                ch_spectre_test_reference,
                params.spectre_test_clair3_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )
            ch_spectre_vcf = call_cnv_spectre.out.vcf
            ch_versions = ch_versions.mix(call_cnv_spectre.out.versions)
        }

        else {

        ch_combined = ch_snv_vcf
        .join(mosdepth.out.regions_bed, by: 0)
        // Result: [meta, vcf_file, bed_file]

        // Transform for cnv_subworkflow - assuming it expects separate channels
        ch_spectre_bed = ch_combined.map { meta, vcf, bed -> bed }
        ch_spectre_vcf = ch_combined.map { meta, vcf, bed -> vcf }

        ch_spectre_reference = ch_samplesheet
        .map { meta, data -> meta.id }
        .join(
            ch_combined.map { meta, vcf, bed -> [meta.id, vcf] },
            by: 0
        )  // Combine with VCF
        .combine(Channel.fromPath(params.fasta_file, checkIfExists: true))
        .map { sample_id, vcf_file, fasta ->
            def meta = [id: sample_id]
            tuple(meta, fasta)
        }

        call_cnv_spectre(
        ch_spectre_bed,
        ch_spectre_reference,
        ch_spectre_vcf,
        params.spectre_metadata,
        params.spectre_blacklist
        )

        ch_spectre_vcf = call_cnv_spectre.out.vcf
        ch_versions = ch_versions.mix(call_cnv_spectre.out.versions)
        }

    }

    if (params.cnv_hificnv){

        call_hificnv(
            ch_input_bam,
            ch_fasta,
            params.exclude_bed_hificnv
        )
        ch_hificnv_vcf = call_hificnv.out.vcf
        ch_versions = ch_versions.mix(call_hificnv.out.versions)
    }

    ch_cnv_vcf = params.cnv_spectre ? ch_spectre_vcf :
            params.cnv_hificnv ? ch_hificnv_vcf :
            Channel.empty()



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

        annotate_str(
            call_str.out.vcf,
            ch_variant_catalogue
        )

        ch_str_vcf = annotate_str.out.vcf
        ch_versions = ch_versions.mix(call_str.out.versions)
    } else {
        ch_str_vcf = Channel.empty()
    }

/*
================================================================================
                            VCF UNIFICATION
================================================================================
*/

if (params.unify_geneyx) {
    ch_sv_unify = ch_sv_vcf_final.map { meta, path ->
    [[id:meta.id], path]  // Extract just the ID string and keep the path
    }

    ch_combined = ch_sv_unify
        .join(ch_cnv_vcf, by: 0, remainder: true)
        .join(ch_str_vcf, by: 0, remainder: true)

    unify_vcf_subworkflow(
        ch_combined.map { meta, sv, cnv, str -> [meta, sv ?: []] },           // ch_sv_vcfs
        ch_combined.map { meta, sv, cnv, str -> [meta, cnv ?: []] },    // ch_cnv_vcf
        ch_combined.map { meta, sv, cnv, str -> [meta, str ?: []] },    // ch_repeat_vcf
        params.modify_str_calls ?: false                                // modify_repeats
    )

    ch_versions = ch_versions.mix(unify_vcf_subworkflow.out.versions)
}
    softwareVersionsToYAML(ch_versions)
    .collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_longraredisease_software_versions.yml',
        sort: true,
        newLine: true
    ).set { ch_collated_versions }

    // ch_collated_versions = Channel.empty()
emit:
    versions = ch_collated_versions


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
