include { LONGPHASE_PHASE } from '../../modules/nf-core/longphase/phase/main.nf'


workflow longphase_variants {
    take:
    ch_bam // channel: tuple(val(meta), path(bam), path(bai))
    ch_snv_vcf
    ch_sv_vcf
    ch_fasta // channel: tuple(val(meta2), path(fasta))
    ch_fai   // channel: tuple(val(meta3), path(fai))

    main:
    ch_versions = Channel.empty()

    ch_input = ch_bam
        .join(ch_snv_vcf, by: 0)
        .join(ch_sv_vcf, by: 0, remainder: true)
        .map { meta, bam, bai, snp_vcf, sv_vcf ->
            // Handle case where sv_vcf might be null
            def sv = sv_vcf ?: []
            tuple(meta, bam, bai, snp_vcf, sv, [])
        }


    LONGPHASE_PHASE(
        ch_input,  // tuple(meta, bam, bai)
        ch_fasta, // tuple(meta2, fasta)
        ch_fai    // tuple(meta3, fai)
    )

    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)

    emit:
    snv_vcf = LONGPHASE_PHASE.out.snv_vcf
    sv_vcf= LONGPHASE_PHASE.out.sv_vcf
    versions = ch_versions
}
