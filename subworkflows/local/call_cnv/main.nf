/*
========================================================================================
    CALL CNV SUBWORKFLOW
========================================================================================
*/

include { SPECTRE_CNVCALLER                   } from '../../../modules/local/spectre/cnvcaller/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SPECTRE } from '../../../modules/nf-core/bcftools/sort'
include { HIFICNV                             } from '../../../modules/nf-core/hificnv/main'

workflow CALL_CNV {
    take:
    ch_input_bam            // channel: [meta, bam, bai]
    ch_mosdepth_summary     // channel: [meta, summary.txt] - can be empty
    ch_mosdepth_regions_bed // channel: [meta, regions.bed.gz] - can be empty
    ch_mosdepth_regions_csi // channel: [meta, regions.bed.gz.csi] - can be empty
    ch_snv_vcf              // channel: [meta, vcf]
    ch_snv_phased_vcf       // channel: [meta, vcf] - for HiFiCNV MAF
    ch_fasta                // channel: [meta, fasta, fai]

    main:
    ch_versions = Channel.empty()
    ch_cnv_vcf = Channel.empty()
    ch_cnv_tbi = Channel.empty()
    ch_cnv_bed = Channel.empty()
    ch_cnv_spc = Channel.empty()
    ch_cnv_winstats = Channel.empty()
    ch_cnv_txt = Channel.empty()
    ch_hificnv_copynum = Channel.empty()
    ch_hificnv_depth = Channel.empty()
    ch_hificnv_maf = Channel.empty()

    // ===================================================================================
    // PacBio / HiFi CNV Calling
    // ===================================================================================
    if ((params.sequencing_platform == 'pacbio' || params.sequencing_platform == 'hifi') || params.filter_targets) {

        // Prepare BAM + BAI + MAF channel
        ch_bam_bai_maf = ch_input_bam
            .join(ch_snv_phased_vcf.map { meta, vcf -> [[id: meta.id], vcf] }, by: 0)
            .map { meta, bam, bai, maf -> [meta, bam, bai, maf] }

        // Create exclude bed channel (optional)
        ch_exclude = params.hificnv_exclude_bed
            ? Channel.of([[id: 'exclude'], file(params.hificnv_exclude_bed, checkIfExists: true)]).first()
            : Channel.of([[id: 'exclude'], []]).first()

        // Create expected CN bed channel (optional)
        ch_expected_cn = params.hificnv_expected_cn_bed
            ? Channel.of([[id: 'expected_cn'], file(params.hificnv_expected_cn_bed, checkIfExists: true)]).first()
            : Channel.of([[id: 'expected_cn'], []]).first()

        HIFICNV(
            ch_bam_bai_maf,
            ch_fasta,
            ch_exclude,
            ch_expected_cn
        )

        ch_cnv_vcf = HIFICNV.out.vcf
        ch_hificnv_copynum = HIFICNV.out.copynum
        ch_hificnv_depth = HIFICNV.out.depth
        ch_hificnv_maf = HIFICNV.out.maf
        ch_versions = ch_versions.mix(HIFICNV.out.versions)
    }

    // ===================================================================================
    // ONT CNV Calling (Spectre)
    // ===================================================================================
    else if (params.sequencing_platform == 'ont' && !params.filter_targets) {

        if (params.use_test_data) {
            // Use test data
            ch_test_meta = Channel.of([id: 'test'])
            ch_test_summary = ch_test_meta.map { meta -> [meta, file(params.spectre_test_summary_txt)] }
            ch_test_regions_bed = ch_test_meta.map { meta -> [meta, file(params.spectre_test_regions_bed)] }
            ch_test_regions_csi = ch_test_meta.map { meta -> [meta, file(params.spectre_test_regions_csi)] }
            ch_test_vcf = ch_test_meta.map { meta -> [meta, file(params.spectre_test_clair3_vcf)] }
            ch_test_fasta = ch_test_meta.map { meta -> [meta, file(params.spectre_test_fasta_file)] }

            ch_mosdepth_dir = ch_test_summary
                .join(ch_test_regions_bed)
                .join(ch_test_regions_csi)
                .map { meta, summary, bed, csi -> [meta, summary, bed, csi] }

            ch_spectre_input = ch_mosdepth_dir
                .join(ch_test_vcf)
                .map { meta, summary, bed, csi, vcf -> [meta, summary, bed, csi, vcf] }

            SPECTRE_CNVCALLER(
                ch_spectre_input,
                ch_test_fasta,
                params.spectre_metadata,
                params.spectre_blacklist,
                params.spectre_bin_size ?: 1000
            )
        } else {
            // Use real mosdepth output
            ch_mosdepth_dir = ch_mosdepth_summary
                .join(ch_mosdepth_regions_bed)
                .join(ch_mosdepth_regions_csi)
                .map { meta, summary, bed, csi -> [meta, summary, bed, csi] }

            ch_spectre_input = ch_mosdepth_dir
                .join(ch_snv_vcf)
                .map { meta, summary, bed, csi, vcf -> [meta, summary, bed, csi, vcf] }

            SPECTRE_CNVCALLER(
                ch_spectre_input,
                ch_fasta,
                params.spectre_metadata,
                params.spectre_blacklist,
                params.spectre_bin_size ?: 1000
            )
        }

        // ch_versions = ch_versions.mix(SPECTRE_CNVCALLER.out.versions)
        BCFTOOLS_SORT_SPECTRE(SPECTRE_CNVCALLER.out.vcf)

        ch_cnv_vcf = BCFTOOLS_SORT_SPECTRE.out.vcf
        ch_cnv_tbi = BCFTOOLS_SORT_SPECTRE.out.tbi
        ch_cnv_bed = SPECTRE_CNVCALLER.out.bed
        ch_cnv_spc = SPECTRE_CNVCALLER.out.spc
        ch_cnv_winstats = SPECTRE_CNVCALLER.out.winstats
        ch_cnv_txt = SPECTRE_CNVCALLER.out.txt
    }

    emit:
    vcf             = ch_cnv_vcf
    tbi             = ch_cnv_tbi
    bed             = ch_cnv_bed
    spc             = ch_cnv_spc
    winstats        = ch_cnv_winstats
    txt             = ch_cnv_txt
    hificnv_copynum = ch_hificnv_copynum
    hificnv_depth   = ch_hificnv_depth
    hificnv_maf     = ch_hificnv_maf
    versions        = ch_versions
}
