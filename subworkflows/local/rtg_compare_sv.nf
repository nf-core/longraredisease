include { RTG_MENDELIAN  as RTG_MENDELIAN_SV              } from '../../modules/local/rtg/mendelian/main.nf'
include { RTG_NONMENDELIAN    as RTG_NONMENDELIAN_SV               } from '../../modules/local/rtg/mendelian_violations/main.nf'

workflow rtg_compare_sv {

    take:
    ch_sdf
    ch_vcf
    ch_ped



    main:
    ch_versions = Channel.empty()


        RTG_MENDELIAN_SV(
            ch_vcf,
            ch_sdf,
            ch_ped,

        )
        ch_mendelian_vcf = RTG_MENDELIAN_SV.out.vcf
        ch_versions = ch_versions.mix(RTG_MENDELIAN_SV.out.versions)



        RTG_NONMENDELIAN_SV(
            ch_vcf,
            ch_sdf,
            ch_ped,

        )

        ch_denovo_vcf = RTG_NONMENDELIAN_SV.out.vcf
        ch_versions = ch_versions.mix(RTG_NONMENDELIAN_SV.out.versions)



    emit:
    mendelian_vcf_sv = ch_mendelian_vcf
    denovo_vcf_sv = ch_denovo_vcf
    versions = ch_versions
}
