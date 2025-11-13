include { HIFICNV                                 } from '../../modules/local/hificnv/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_HIFICNV } from '../../modules/nf-core/bcftools/sort'

workflow call_hificnv {
    take:
    ch_mosdepth_output    // channel: [ val(meta), path(mosdepth_dir) ]
    ch_reference          // channel: [ val(meta2), path(fasta) ]
    ch_exclude_bed        // path to blacklist file

    main:

    //
    // MODULE: Run HiFiCNV CNV calling
    //
    HIFICNV(
        ch_mosdepth_output,
        ch_reference,
        ch_exclude_bed
    )

    BCFTOOLS_SORT_HIFICNV(HIFICNV.out.vcf)

    emit:
    vcf       = BCFTOOLS_SORT_HIFICNV.out.vcf
    tbi       = BCFTOOLS_SORT_HIFICNV.out.tbi
    bedgraph  = HIFICNV.out.cnval
    versions  = HIFICNV.out.versions

}
