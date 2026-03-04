include { RTGTOOLS_FORMAT                   } from '../../modules/local/rtgtools/format_reference/main.nf'
include { RTGTOOLS_MENDELIAN                } from '../../modules/local/rtgtools/vcf_comparison_mendelian/main.nf'
include { RTGTOOLS_DENOVO                   } from '../../modules/local/rtgtools/vcf_comparison_denovo/main.nf'

workflow rtg_trio_comparison_subworkflow {

    take:
    ch_fasta
    ch_vcf
    ch_ped
    run_mendelian
    run_denovo

    main:
    ch_versions = Channel.empty()

    RTGTOOLS_FORMAT(ch_fasta)

    // Conditionally run based on params
    if (run_mendelian) {
        RTGTOOLS_MENDELIAN(
            ch_vcf,
            RTGTOOLS_FORMAT.out.sdf,
            ch_ped,

        )
        ch_mendelian_vcf = RTGTOOLS_MENDELIAN.out.vcf
    } else {
        ch_mendelian_vcf = Channel.empty()
    }

    if (run_denovo) {
        RTGTOOLS_DENOVO(
            ch_vcf,
            RTGTOOLS_FORMAT.out.sdf,
            ch_ped,

        )
        ch_denovo_vcf = RTGTOOLS_DENOVO.out.vcf
    } else {
        ch_denovo_vcf = Channel.empty()
    }


    emit:
    sdf = RTGTOOLS_FORMAT.out.sdf
    mendelian_vcf = ch_mendelian_vcf
    denovo_vcf = ch_denovo_vcf
}
