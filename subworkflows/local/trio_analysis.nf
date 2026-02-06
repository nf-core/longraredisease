include { CREATE_PEDIGREE_FILE } from '../../modules/local/create_ped_file/main.nf'
include { GLNEXUS_TRIO } from '../../modules/local/glnexus_trio/main.nf'
include { SNIFFLES2_TRIO } from '../../modules/local/sniffles2_trio/main.nf'
include { PHASE_TRIO } from '../../modules/local/whatshap_phase_trio/main.nf'
include { RTG_MENDELIAN } from '../../modules/local/rtg_mendelian/main.nf'

workflow trio_analysis {
    take:
    ch_trio_bams        // channel: [family_meta, [proband_bam, pat_bam, mat_bam], [proband_bai, pat_bai, mat_bai]]
    ch_fasta           // channel: [meta, fasta]
    ch_fai             // channel: [meta, fai]
    ch_pedigree        // channel: [family_meta, pedigree_file]
    ch_trf             // channel: [meta, tandem_repeat_bed]

    main:
    ch_versions = Channel.empty()

    // Create PED from samplesheet
        //
        ch_input
            .map { meta, _files -> [ [ id: meta.project ], meta ] }
            .groupTuple()
            .set { ch_samplesheet_ped_in }

        SAMPLESHEET_PED ( ch_samplesheet_ped_in )
        ch_versions = ch_versions.mix(SAMPLESHEET_PED.out.versions)

        SAMPLESHEET_PED.out.ped
            .collect()
            .set { ch_samplesheet_pedfile }

    // Small variant calling for trio
    if (params.snv) {
        CLAIR3_NOVA_TRIO(
            ch_trio_bams,
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(CLAIR3_NOVA_TRIO.out.versions)

        // Joint genotyping
        GLNEXUS_TRIO(
            CLAIR3_NOVA_TRIO.out.gvcf,
            ch_pedigree
        )
        ch_versions = ch_versions.mix(GLNEXUS_TRIO.out.versions)

        ch_trio_snv_vcf = GLNEXUS_TRIO.out.vcf
    } else {
        ch_trio_snv_vcf = Channel.empty()
    }

    // Structural variant calling for trio
    if (params.sv) {
        SNIFFLES2_TRIO(
            ch_trio_bams,
            ch_fasta,
            ch_trf
        )
        ch_versions = ch_versions.mix(SNIFFLES2_TRIO.out.versions)

        ch_trio_sv_vcf = SNIFFLES2_TRIO.out.vcf
    } else {
        ch_trio_sv_vcf = Channel.empty()
    }

    // Pedigree phasing
    if (params.phase && params.snv) {
        WHATSHAP_PHASE_TRIO(
            ch_trio_snv_vcf,
            ch_trio_bams,
            ch_pedigree
        )
        ch_versions = ch_versions.mix(WHATSHAP_PHASE_TRIO.out.versions)

        ch_phased_vcf = WHATSHAP_PHASE_TRIO.out.vcf
        ch_haplotagged_bams = WHATSHAP_PHASE_TRIO.out.haplotagged_bams
    } else {
        ch_phased_vcf = ch_trio_snv_vcf
        ch_haplotagged_bams = Channel.empty()
    }

    // Mendelian inheritance analysis
    if (params.snv || params.sv) {
        RTG_MENDELIAN(
            ch_trio_snv_vcf.mix(ch_trio_sv_vcf),
            ch_pedigree
        )
        ch_versions = ch_versions.mix(RTG_MENDELIAN.out.versions)
    }

    emit:
    snv_vcf = ch_trio_snv_vcf
    sv_vcf = ch_trio_sv_vcf
    phased_vcf = ch_phased_vcf
    haplotagged_bams = ch_haplotagged_bams
    mendelian_report = RTG_MENDELIAN.out.report
    versions = ch_versions
}
