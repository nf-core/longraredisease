// Copy Number Variant Detection Subworkflow

include { SPECTRE_CNVCALLER } from '../../modules/local/spectre/cnvcaller/main'
include { SPECTRE_ROUNDDP } from '../../modules/local/spectre/rounddp/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SPECTRE } from '../../modules/nf-core/bcftools/sort'



workflow call_spectre_cnv {
    take:
    ch_mosdepth_output    // channel: [ val(meta), path(mosdepth_dir) ]
    ch_reference          // channel: [ val(meta2), path(fasta) ]
    ch_snv_vcf            // channel: [ val(meta3), path(vcf) ]
    ch_metadata          // path to metadata file
    ch_blacklist         // path to blacklist file

    main:

    //
    // MODULE: Run SPECTRE CNV calling
    //
    ch_versions = Channel.empty()

    SPECTRE_CNVCALLER(
        ch_mosdepth_output,
        ch_reference,
        ch_snv_vcf,
        ch_metadata,
        ch_blacklist
    )

    ch_versions = ch_versions.mix(SPECTRE_CNVCALLER.out.versions)

    SPECTRE_ROUNDDP(SPECTRE_CNVCALLER.out.vcf)

    ch_versions = ch_versions.mix(SPECTRE_ROUNDDP.out.versions)

    BCFTOOLS_SORT_SPECTRE(SPECTRE_ROUNDDP.out.vcf)

    ch_versions = ch_versions.mix(BCFTOOLS_SORT_SPECTRE.out.versions)

    emit:
    vcf       = BCFTOOLS_SORT_SPECTRE.out.vcf
    tbi       = BCFTOOLS_SORT_SPECTRE.out.tbi
    bed       = SPECTRE_CNVCALLER.out.bed
    bed_index = SPECTRE_CNVCALLER.out.bed_index
    spc       = SPECTRE_CNVCALLER.out.spc
    winstats  = SPECTRE_CNVCALLER.out.winstats
    versions  = ch_versions

}
