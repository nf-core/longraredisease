include { LONGPHASE_HAPLOTAG } from '../../../modules/nf-core/longphase/haplotag/main.nf'


workflow HAPLOTAG_BAM {
    take:
    ch_bam // channel: tuple(val(meta), path(bam), path(bai))
    ch_phased_snv
    ch_phased_sv
    ch_fasta // channel: tuple(val(meta2), path(fasta))
    ch_fai   // channel: tuple(val(meta3), path(fai))

    main:
    ch_versions = Channel.empty()

    ch_input = ch_bam
        .join(ch_phased_snv, by: 0)
        .join(ch_phased_sv, by: 0, remainder: true)
        .map { meta, bam, bai, snp_vcf, sv_vcf ->
            // Handle case where sv_vcf might be null
            def sv = sv_vcf ?: []
            tuple(meta, bam, bai, snp_vcf, sv, [])
        }

    LONGPHASE_HAPLOTAG(
        ch_input,  // tuple(meta, bam, bai)
        ch_fasta, // tuple(meta2, fasta)
        ch_fai    // tuple(meta3, fai)
    )

    ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

    emit:
    bam = LONGPHASE_HAPLOTAG.out.bam
    versions = ch_versions
}
