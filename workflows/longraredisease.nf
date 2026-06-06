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
include { BAM2FASTQ                          } from '../subworkflows/local/bam2fastq/main.nf'
include { ALIGN                              } from '../subworkflows/local/align/main.nf'
include { CAT_FASTQ                          } from '../modules/nf-core/cat/fastq/main.nf'
include { NANOPLOT as NANOPLOT_QC            } from '../modules/nf-core/nanoplot/main'
include { CREATE_PEDIGREE_FILE               } from '../modules/local/create_ped_file/main.nf'

// Coverage analysis subworkflows
include { MOSDEPTH_SUBWORKFLOW               } from '../subworkflows/local/mosdepth/main.nf'
include { MULTIQC_MOSDEPTH                   } from '../modules/local/multiqc_mosdepth/main.nf'
include { MULTIQC                            } from '../modules/nf-core/multiqc/main.nf'

// Trio analysis - rtg format reference file
include { RTGTOOLS_FORMAT                    } from '../modules/nf-core/rtgtools/format/main.nf'

// Methylation calling
include { METHYL                             } from '../subworkflows/local/methyl/main.nf'

// SNV/indel calling
include { CALL_SNV                           } from '../subworkflows/local/call_snv/main.nf'
include { ANNOTATE_SNV                       } from '../subworkflows/local/annotate_snv/main.nf'

// Haplotag BAM
include { SNIFFLES as SNIFFLES_UNPHASED      } from '../modules/nf-core/sniffles/main.nf'
include { LONGPHASE_PHASE                    } from '../modules/nf-core/longphase/phase/main.nf'
include { HAPLOTAG_BAM                       } from '../subworkflows/local/haplotag_bam/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HAPLOTAG } from '../modules/nf-core/samtools/index/main'

// SV calling
include { CALL_SV                            } from '../subworkflows/local/call_sv/main.nf'
include { SNIFFLES_GENERATE_PLOTS            } from '../modules/local/sniffles/generate_plots/main.nf'
include { FILTER_SV as FILTER_SV_SNIFFLES    } from '../subworkflows/local/filter_sv/main.nf'

// Annotate and prioritize variants
include { SVANNA_PRIORITIZE                  } from '../modules/local/svanna/main.nf'

// SV calling for trios
include { TRIO_CONCORDANCE_SV                     } from '../subworkflows/local/trio_concordance_sv/main.nf'


// SNV calling for trios
include { JOINT_GENOTYPE_SNV                 } from '../subworkflows/local/joint_genotype_snv/main.nf'
include { TRIO_CONCORDANCE_SNV                    } from '../subworkflows/local/trio_concordance_snv/main.nf'

// STR analysis subworkflow
include { CALL_STR                          } from '../subworkflows/local/call_str/main.nf'

// CNV calling subworkflows
include { CALL_CNV                           } from '../subworkflows/local/call_cnv/main.nf'

// Merge SV - multiple callers
include { FILTER_SV  as FILTER_SV_SVIM       } from '../subworkflows/local/filter_sv/main.nf'
include { ANNOTATE_SV                        } from '../subworkflows/local/annotate_sv/main.nf'
include { FILTER_SV  as FILTER_SV_CUTESV     } from '../subworkflows/local/filter_sv/main.nf'
include { GUNZIP as GUNZIP_SVIM              } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_CUTESV            } from '../modules/nf-core/gunzip/main.nf'
include { MERGE_SV                           } from '../subworkflows/local/merge_sv/main.nf'

include { UNIFYVCF                          } from '../modules/local/unify_vcf/main.nf'
include { ANNOTATE_UNIFIED                   } from '../subworkflows/local/annotate_unified/main.nf'

// VCF processing subworkflows
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_longraredisease_pipeline'
include { paramsSummaryMap                   } from 'plugin/nf-schema'
include { paramsSummaryMultiqc               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { citationBibliographyText           } from '../subworkflows/local/utils_nfcore_longraredisease_pipeline'


workflow LONGRAREDISEASE {

    main:

    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .map { meta, file_path ->
            [meta, file(file_path)]
        }

    if (params.trio_analysis) {
        pedigree_input = ch_samplesheet
            .map { meta, file_path ->
                [meta.family_id, meta, file_path]
            }
            .groupTuple(by: 0, size: 3)
            .map { family_id, metas, file_paths ->
                def family_meta = [id: family_id]
                [family_meta, metas, file_paths]
            }

        CREATE_PEDIGREE_FILE(
            pedigree_input.map { family_meta, metas, file_paths ->
                [family_meta, metas]
            }
        )
    }

/*
=======================================================================================
                                REFERENCE FILES SETUP
=======================================================================================
*/
    // Initialize versions channel
    ch_versions = channel.empty()

    def fasta_path = params.genome
        ? params.genomes.containsKey(params.genome)
            ? params.genomes[params.genome].fasta
            : error("Genome '${params.genome}' not found in iGenomes config. Check --genome or set --igenomes_ignore false.")
        : params.fasta_file

    ch_fasta = Channel
        .fromPath(fasta_path, checkIfExists: true)
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

        RTGTOOLS_FORMAT(ch_fasta.map { meta, fasta -> [meta, fasta, [], []] })

        // Extract unique family IDs from samplesheet
        ch_family_ids = ch_samplesheet
            .map { meta, file_path -> meta.family_id }
            .unique()

        // Replicate SDF with family-specific metadata
        ch_sdf = RTGTOOLS_FORMAT.out.sdf
            .map { meta, sdf -> sdf }
            .combine(ch_family_ids)
            .map { sdf, family_id -> [[id: family_id], sdf] }
    }

    // Tandem repeat file for Sniffles (only if SV calling is enabled)
    if (params.sv) {
        ch_trf = Channel
            .fromPath(params.sniffles_tandem_file, checkIfExists: true)
            .map { bed -> tuple([id: "trf"], bed) }
            .first()
    } else {
        ch_trf = channel.empty()
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
            .map { meta, file_path ->
            // file_path is already a Path object from the ch_samplesheet .map { meta, fp -> [meta, file(fp)] }
            if (file_path.isFile() && (file_path.name.endsWith('.fastq.gz') || file_path.name.endsWith('.fq.gz'))) {
                // Single FASTQ file case
                return [meta, [file_path]]
            } else if (file_path.isDirectory()) {
                // Directory with multiple FASTQ files
                def fastq_files = file_path.listFiles().findAll {
                    it.name.endsWith('.fastq.gz') || it.name.endsWith('.fq.gz')
                }

                if (fastq_files.isEmpty()) {
                    error "No FASTQ files found in directory: ${file_path} for sample ${meta.id}"
                }

                return [meta, fastq_files]
            } else {
                error "Invalid FASTQ input for sample ${meta.id}: ${file_path}"
            }
        }

        ch_fastq_files.branch { meta, files ->
            single: files.size() == 1
                return [meta, files[0]]
            multiple: files.size() > 1
                return [meta + [single_end: true], files]
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
        ALIGN (
            ch_fasta,
            ch_processed_fastq,
            params.winnowmap_kmers,
            params.filter_targets
        )

        ch_versions = ch_versions.mix(ALIGN.out.versions)

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = ALIGN.out.bam
        ch_final_sorted_bai = ALIGN.out.bai

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
            .map { meta, file_path ->
                if (!file_path) {
                    error "No BAM input provided for sample ${meta.id}"
                }

                if (file_path.isFile() && file_path.name.endsWith('.bam')) {
                    // Single BAM file case
                    return [meta + [is_multiple: false], file_path]
                } else if (file_path.isDirectory()) {
                    // Directory with multiple BAM files case
                    def bam_files = file_path.listFiles().findAll { it.name.endsWith('.bam') }

                    if (bam_files.isEmpty()) {
                        error "No BAM files found for sample ${meta.id} in directory: ${file_path}"
                    }

                    return [meta + [is_multiple: bam_files.size() > 1], bam_files]
                } else {
                    error "Invalid BAM input for sample ${meta.id}: ${file_path} (not a file or directory)"
                }
            }

        // Convert BAM to FASTQ
        BAM2FASTQ (
            ch_bam_files,
            [[:], []],
            [[:], []],
            [[:], []]
        )

        ch_versions = ch_versions.mix(BAM2FASTQ.out.versions)
        // Align FASTQ reads to reference genome using minimap2
        ALIGN (
            ch_fasta,
            BAM2FASTQ.out.other,
            params.winnowmap_kmers,
            params.filter_targets
        )

        ch_versions = ch_versions.mix(ALIGN.out.versions)

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = ALIGN.out.bam
            .map { meta, bam ->
                def clean_meta = [id: meta.id]
                [clean_meta, bam]
            }

        ch_final_sorted_bai = ALIGN.out.bai
            .map { meta, bai ->
                def clean_meta = [id: meta.id]
                [clean_meta, bai]
            }


        // Prepare input for nanoplot from FASTQ
        ch_nanoplot = BAM2FASTQ.out.other
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
            .map { meta, file_path ->
                // file_path is already a Path object (file() was applied earlier in ch_samplesheet)
                def bai_file = file("${file_path}.bai", checkIfExists: true)
                return [meta, file_path, bai_file]
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
            tuple(meta, bam, bai, [])
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

        METHYL(
            ch_methyl_input,
            ch_fasta_fai,
            [[:], []]
        )

        ch_versions = ch_versions.mix(METHYL.out.versions)
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

        MOSDEPTH_SUBWORKFLOW(
            ch_input_bam_bai_bed,
            [[:], []]
        )
        ch_versions = ch_versions.mix(MOSDEPTH_SUBWORKFLOW.out.versions)

        // Combine all mosdepth outputs per sample, preserving metadata

        ch_mosdepth = MOSDEPTH_SUBWORKFLOW.out.global_txt
            .join(MOSDEPTH_SUBWORKFLOW.out.summary_txt)
            .join(MOSDEPTH_SUBWORKFLOW.out.regions_txt)
            .map { meta, file1, file2, file3 ->
                [meta, [file1, file2, file3]]
            }

        MULTIQC_MOSDEPTH (
            ch_mosdepth  // Pass [meta, [files]] tuples
        )

        ch_versions = ch_versions.mix(MULTIQC_MOSDEPTH.out.versions)
    }


/*
=======================================================================================
                                CALL SNV/INDEL
=======================================================================================
*/

    if (params.snv || params.haplotag_bam) {

        CALL_SNV (
            ch_input_bam,
            ch_fasta,
            ch_fai,
            params.run_deepvariant,
            ch_input_bam_bai_bed
        )

        ch_versions = ch_versions.mix(CALL_SNV.out.versions)

        ch_snv_vcf = CALL_SNV.out.vcf
        ch_snv_phased_vcf = CALL_SNV.out.phased_vcf
    }

    if (params.snv && params.annotate_clair3) {
        ANNOTATE_SNV(
            ch_snv_vcf,
            params.snpeff_db
        )
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

        ch_input_longphase = ch_input_bam  // [meta, bam, bai]
        .join(ch_snv_vcf, by: 0)
        .join(SNIFFLES_UNPHASED.out.vcf, by: 0, remainder: true)
        .map { meta, bam, bai, snv_vcf, sv_vcf ->
             def sv = sv_vcf ?: []
             tuple(meta, bam, bai, snv_vcf, sv, [])  // [] = no mod file
        }

        LONGPHASE_PHASE(
            ch_input_longphase,
            ch_fasta,
            ch_fai
        )

        HAPLOTAG_BAM(
            ch_input_bam,
            LONGPHASE_PHASE.out.snv_vcf,
            LONGPHASE_PHASE.out.sv_vcf,
            ch_fasta,
            ch_fai
        )

        SAMTOOLS_INDEX_HAPLOTAG(HAPLOTAG_BAM.out.bam)

        ch_input_bam = HAPLOTAG_BAM.out.bam
        .join(SAMTOOLS_INDEX_HAPLOTAG.out.bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }

        ch_versions = ch_versions.mix(SNIFFLES_UNPHASED.out.versions)
        ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)
        ch_versions = ch_versions.mix(HAPLOTAG_BAM.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HAPLOTAG.out.versions)
    }


/*
=======================================================================================
                                CALL SV
=======================================================================================
*/
    if (params.sv){

        CALL_SV(
            ch_input_bam,
            ch_fasta,
            ch_trf,
            params.vcf_output,
            params.snf_output,
            params.merge_sv,
            params.run_svim
        )

        ch_sv_vcf_final = CALL_SV.out.sniffles_vcf
        ch_svim_vcf = CALL_SV.out.svim_vcf
        ch_versions = ch_versions.mix(CALL_SV.out.versions)

        if (params.filter_pass_sv) {

            FILTER_SV_SNIFFLES(
                CALL_SV.out.sniffles_vcf_tbi
                    .filter { _meta, vcf, _tbi -> vcf != null }
                    .map { meta, vcf, tbi -> [meta + [caller: 'sniffles'], vcf, tbi] },
                params.coverage_bed,
                params.downsample_sv,
                MOSDEPTH_SUBWORKFLOW.out.summary_txt,
                MOSDEPTH_SUBWORKFLOW.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support,
                params.min_read_support_limit
            )
            ch_sv_vcf_final = FILTER_SV_SNIFFLES.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
            ch_versions = ch_versions.mix(FILTER_SV_SNIFFLES.out.versions)

        }

        if (params.filter_pass_sv && params.run_svim) {
            FILTER_SV_SVIM(
                CALL_SV.out.svim_vcf_tbi
                    .filter { _meta, vcf, _tbi -> vcf != null }
                    .map { meta, vcf, tbi -> [meta + [caller: 'svim'], vcf, tbi] },
                params.coverage_bed,
                params.downsample_sv,
                MOSDEPTH_SUBWORKFLOW.out.summary_txt,
                MOSDEPTH_SUBWORKFLOW.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support,
                params.min_read_support_limit
            )

            ch_svim_vcf = FILTER_SV_SVIM.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
        }

    }

    else {
        ch_sv_vcf_final = channel.empty()
    }

    /*
    ================================================================================
                            Annotate SV vcf
    ================================================================================
    */

    if (params.sv && params.annotate_sv) {
        ANNOTATE_SV(
            ch_samplesheet,
            ch_sv_vcf_final,
            ch_svim_vcf,
            ch_snv_vcf,
            [],
            [],
            []
        )
    }

/*
================================================================================
                        Sniffles SV ANNOTATION WITH SVANNA
================================================================================
*/

    if (params.sv && params.run_svanna) {
        // Filter samplesheet to only include samples with HPO terms
        ch_samplesheet_with_hpo = ch_samplesheet
            .filter { meta, _file_path ->
                meta.hpo_terms && meta.hpo_terms.trim() != ""
            }

        ch_hpo_terms = ch_samplesheet_with_hpo
            .map { meta, _file_path ->
                [meta.id, meta.hpo_terms]
            }

        ch_sv_svanna = ch_sv_vcf_final
            .map { meta, vcf -> [meta.id, meta, vcf] }
            .join(ch_hpo_terms, by: 0)
            .map { _sample_id, meta, vcf, hpo_terms ->
                def meta_with_hpo = meta + [hpo_terms: hpo_terms]
                [meta_with_hpo, vcf, hpo_terms]
            }


        SVANNA_PRIORITIZE(
            ch_sv_svanna.map { meta, vcf, hpo_terms -> [meta, vcf] },
            params.svanna_db,
            ch_sv_svanna.map { meta, vcf, hpo_terms -> hpo_terms }
        )

        ch_versions = ch_versions.mix(SVANNA_PRIORITIZE.out.versions)
    }
/*
=======================================================================================
                                Trio analysis
=======================================================================================
*/

    if (params.sv && params.trio_analysis) {
        TRIO_CONCORDANCE_SV(
            ch_sdf,
            CALL_SV.out.sniffles_snf,
            ch_samplesheet,
            ch_fasta,
            CREATE_PEDIGREE_FILE.out.ped
                .map { meta, ped -> [meta, ped] }
        )
    }

    if (params.snv && params.trio_analysis) {

        ch_gvcf = CALL_SNV.out.gvcf

        JOINT_GENOTYPE_SNV(
            ch_gvcf,
            ch_samplesheet,
            [[:], []]
        )

        ch_versions = ch_versions.mix(JOINT_GENOTYPE_SNV.out.versions)

        ch_trio_snv_vcf = JOINT_GENOTYPE_SNV.out.vcf
            .map { meta, vcf -> [meta + [variant_type: 'snv'], vcf] }

        TRIO_CONCORDANCE_SNV(
            ch_sdf,
            ch_trio_snv_vcf,
            CREATE_PEDIGREE_FILE.out.ped
                .map { meta, ped -> [meta, ped] }
        )
    }

    // annotate trios - future release


/*
===============================================================================
                        SHORT TANDEM REPEAT ANALYSIS
================================================================================
*/

    if (params.str) {
        CALL_STR(
            ch_input_bam,
            ch_fasta,
            ch_fai
        )

        ch_str_vcf  = CALL_STR.out.vcf
        ch_versions = ch_versions.mix(CALL_STR.out.versions)
    } else {
        ch_str_vcf = channel.empty()
    }

/*
=======================================================================================
                        COPY NUMBER VARIANT CALLING
=======================================================================================
*/

    if (params.cnv) {
        CALL_CNV(
            ch_input_bam,
            params.sequencing_platform == 'ont' && !params.filter_targets ? MOSDEPTH_SUBWORKFLOW.out.summary_txt : channel.empty(),
            params.sequencing_platform == 'ont' && !params.filter_targets ? MOSDEPTH_SUBWORKFLOW.out.regions_bed : channel.empty(),
            params.sequencing_platform == 'ont' && !params.filter_targets ? MOSDEPTH_SUBWORKFLOW.out.regions_csi : channel.empty(),
            ch_snv_vcf,
            ch_snv_phased_vcf,
            ch_fasta
        )

        ch_cnv_vcf = CALL_CNV.out.vcf
        ch_versions = ch_versions.mix(CALL_CNV.out.versions)
    } else {
        ch_cnv_vcf = channel.empty()
    }

/*
================================================================================
                            MERGE SV with Jasmine
================================================================================
*/

    // Gunzip VCFs for Jasmine (requires uncompressed input)
    if (params.sv && params.merge_sv) {
        if (params.filter_pass_sv) {
            FILTER_SV_CUTESV(
                CALL_SV.out.cutesv_vcf_tbi
                    .filter { _meta, vcf, _tbi -> vcf != null }
                    .map { meta, vcf, tbi -> [meta + [caller: 'cutesv'], vcf, tbi] },
                params.coverage_bed,
                params.downsample_sv,
                MOSDEPTH_SUBWORKFLOW.out.summary_txt,
                MOSDEPTH_SUBWORKFLOW.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support,
                params.min_read_support_limit
            )

            ch_cutesv_vcf = FILTER_SV_CUTESV.out.ch_vcf_tbi.map { meta, vcf, tbi -> [meta, vcf] }
        }
        // Jasmine requires unzipped VCFs

        GUNZIP_SVIM(ch_svim_vcf)

        GUNZIP_CUTESV(ch_cutesv_vcf)

        ch_versions = ch_versions.mix(GUNZIP_SVIM.out.versions)
        ch_versions = ch_versions.mix(GUNZIP_CUTESV.out.versions)

        // Prepare input for JASMINESV - group all uncompressed VCFs by sample
        jasmine_input_ch = CALL_SV.out.sniffles_unzipped_vcf
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
        MERGE_SV(
            jasmine_input_ch,
            ch_fasta,
            ch_fai,
            []
        )

        ch_versions = ch_versions.mix(MERGE_SV.out.versions)

        // Set final SV VCF to merged result for the unify vcf if required
        ch_sv_vcf_final = MERGE_SV.out.vcf
    }


/*
================================================================================
                            VCF UNIFICATION
================================================================================
*/

    if (params.unify_vcf && params.sv && params.cnv && params.str) {
        ch_combined = ch_sv_vcf_final
            .join(ch_cnv_vcf, by: 0, remainder: true)
            .join(ch_str_vcf, by: 0, remainder: true)

        UNIFYVCF(
            ch_combined.map { meta, sv, _cnv, _str -> [meta, sv] },
            ch_combined.map { meta, _sv, cnv, _str -> [meta, cnv ?: []] },
            ch_combined.map { meta, _sv, _cnv, str -> [meta, str ?: []] },
            params.modify_str_calls ?: false
        )

        ch_versions = ch_versions.mix(UNIFYVCF.out.versions)

        if (params.annotate_unified_vcf) {
            ANNOTATE_UNIFIED(UNIFYVCF.out.unified_vcf, params.snpeff_db)
        }
    }

/*
=======================================================================================
                                MULTIQC
=======================================================================================
*/
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(
        Channel.of(paramsSummaryMultiqc(paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")))
            .collectFile(name: 'workflow_summary_mqc.yaml')
    )

    if (params.qc) {
        ch_multiqc_files = ch_multiqc_files.mix(
            NANOPLOT_QC.out.txt.collect { meta, txt -> txt }
        )
    }

    MULTIQC(
        ch_multiqc_files.collect().map { files ->
            [
                [id: 'multiqc'],
                files,
                [file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)],
                file("$projectDir/assets/nf-core-longraredisease_logo_light.png", checkIfExists: true),
                [],
                []
            ]
        }
    )

    //
    // Collate and save software versions
    // Supports both legacy ch_versions and nf-core 3.5+ topic channel versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'longraredisease_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report }
}
