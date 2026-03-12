/*
========================================================================================
    CALL SPECTRE CNV SUBWORKFLOW
========================================================================================
*/

include { SPECTRE_CNVCALLER } from '../../modules/local/spectre/cnvcaller/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SPECTRE } from '../../modules/nf-core/bcftools/sort'


workflow call_spectre_cnv {
    take:
    ch_mosdepth_summary     // channel: [meta, summary.txt]
    ch_mosdepth_regions_bed // channel: [meta, regions.bed.gz]
    ch_mosdepth_regions_csi // channel: [meta, regions.bed.gz.csi]
    ch_snv_vcf              // channel: [meta, vcf]
    ch_fasta                // channel: [meta, fasta]
    ch_metadata             // path to metadata file
    ch_blacklist            // path to blacklist file
    bin_size                // val: bin size for CNV calling

    main:
    ch_versions = channel.empty()


    ch_mosdepth_dir = ch_mosdepth_summary
        .join(ch_mosdepth_regions_bed)
        .join(ch_mosdepth_regions_csi)
        .map { tuple ->
            [tuple[0], tuple[1], tuple[2], tuple[3]]
        }

    ch_input = ch_mosdepth_dir
        .join(ch_snv_vcf)
        .map { tuple ->
            // tuple = [meta, summary, regions_bed, regions_csi, vcf]
            [tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]]
        }

    SPECTRE_CNVCALLER(
        ch_input,
        ch_fasta,
        ch_metadata,
        ch_blacklist,
        bin_size
    )

    ch_versions = ch_versions.mix(SPECTRE_CNVCALLER.out.versions)
    BCFTOOLS_SORT_SPECTRE(SPECTRE_CNVCALLER.out.vcf)


    emit:
    vcf       = BCFTOOLS_SORT_SPECTRE.out.vcf
    tbi       = BCFTOOLS_SORT_SPECTRE.out.tbi
    bed       = SPECTRE_CNVCALLER.out.bed
    spc       = SPECTRE_CNVCALLER.out.spc
    winstats  = SPECTRE_CNVCALLER.out.winstats
    txt       = SPECTRE_CNVCALLER.out.txt
    versions  = ch_versions
}
