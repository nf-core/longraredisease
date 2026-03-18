// Subworkflow to identify STR repeats in a genome using STAGLR

include { STRAGLR } from '../../../modules/local/straglr/main'
include { STRANGER   } from '../../../modules/nf-core/stranger/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_STRAGLR } from '../../../modules/nf-core/bcftools/sort/main.nf'
include { TRGT_GENOTYPE   } from '../../../modules/nf-core/trgt/genotype/main.nf'

workflow CALL_STR {
    take:
    ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta  // channel: [ val(meta2), path(reference) ]
    ch_fai


    main:
    ch_versions = Channel.empty()

    if (params.sequencing_platform== 'ont'){
    STRAGLR(
        ch_bam_bai,
        ch_fasta,
        params.straglr_bed
    )

    BCFTOOLS_SORT_STRAGLR(
        STRAGLR.out.vcf
    )

    ch_variant_catalogue = channel.fromPath(params.variant_catalogue)
        .map { file -> [ [id: 'variant_catalog'], file ] }
        .first()


    STRANGER(BCFTOOLS_SORT_STRAGLR.out.vcf,
    ch_variant_catalogue)

    str_vcf = STRANGER.out.vcf

    ch_versions = ch_versions.mix(STRAGLR.out.versions)}

    if (params.sequencing_platform == 'pacbio'){

        ch_bam_bai_karyo = ch_bam_bai.map {meta, bam, bai -> tuple(meta, bam, bai, [])}

        ch_trgt_bed = channel.fromPath(params.trgt_bed)
        .map { file -> [ [id: 'trgt_bed'], file ] }
        .first()

        TRGT_GENOTYPE(ch_bam_bai_karyo, ch_fasta, ch_fai, ch_trgt_bed )
        str_vcf = TRGT_GENOTYPE.out.vcf
    }


    emit:
    vcf      = str_vcf      // channel: [ val(meta), path(vcf) ]
    versions = ch_versions         // channel: path(versions.yml)
}
