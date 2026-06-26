include { RTG_MENDELIAN  as RTG_MENDELIAN_SNV              } from '../../../modules/local/rtg/mendelian/main.nf'
include { RTG_VIOLATIONS    as RTG_VIOLATIONS_SNV               } from '../../../modules/local/rtg/violations/main.nf'

workflow TRIO_CONCORDANCE_SNV {

    take:
    ch_sdf
    ch_vcf
    ch_ped


    main:
    ch_versions = Channel.empty()

    // Conditionally run based on params
        RTG_MENDELIAN_SNV(
            ch_vcf,
            ch_sdf,
            ch_ped,

        )
        ch_mendelian_vcf = RTG_MENDELIAN_SNV.out.vcf
        ch_versions = ch_versions.mix(RTG_MENDELIAN_SNV.out.versions)



        RTG_VIOLATIONS_SNV(
            ch_vcf,
            ch_sdf,
            ch_ped,

        )
        ch_violations_vcf = RTG_VIOLATIONS_SNV.out.vcf
        ch_versions = ch_versions.mix(RTG_VIOLATIONS_SNV.out.versions)


    emit:
    mendelian_vcf_sv = ch_mendelian_vcf
    denovo_vcf_sv = ch_violations_vcf
    versions = ch_versions
}
