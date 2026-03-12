include { HIFICNV } from '../../modules/nf-core/hificnv/main'

workflow call_hificnv {

    take:
    bam_bai_maf       // channel: tuple val(meta), path(bam), path(bai), path(maf)
    ref_fasta         // channel: tuple val(meta), path(fasta), path(fai)
    exclude_bed       // channel: tuple val(meta), path(bed) - optional exclude regions
    expected_cn_bed   // channel: tuple val(meta), path(bed) - optional expected copy number

    main:
    ch_versions = Channel.empty()

    //
    // Run HiFiCNV with MAF file
    //
    HIFICNV (
        bam_bai_maf,
        ref_fasta,
        exclude_bed,
        expected_cn_bed
    )
    ch_versions = ch_versions.mix(HIFICNV.out.versions)

    emit:
    copynum  = HIFICNV.out.copynum   // channel: tuple val(meta), path(*.copynum.bedgraph)
    depth    = HIFICNV.out.depth     // channel: tuple val(meta), path(*.depth.bw)
    maf      = HIFICNV.out.maf       // channel: tuple val(meta), path(*.maf.bw)
    vcf      = HIFICNV.out.vcf       // channel: tuple val(meta), path(*.vcf.gz)
    versions = ch_versions           // channel: path(versions.yml)
}

