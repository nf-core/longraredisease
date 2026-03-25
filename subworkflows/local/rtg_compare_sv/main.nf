include { SNIFFLES_TRIO                      } from '../../../modules/local/sniffles/trio/main.nf'
include { RTG_MENDELIAN  as RTG_MENDELIAN_SV              } from '../../../modules/local/rtg/mendelian/main.nf'
include { RTG_NONMENDELIAN    as RTG_NONMENDELIAN_SV               } from '../../../modules/local/rtg/mendelian_violations/main.nf'

workflow RTG_COMPARE_SV {

    take:
    ch_sdf
    snf_files
    samplesheet
    fasta
    ch_ped



    main:
    ch_versions = Channel.empty()

    ch_snf_for_trio = snf_files
        .map { meta, snf -> [meta.id, snf] }
        .join(
            samplesheet.map { meta, data -> [data.id, data.family_id] },
            by: 0
        )
        .filter { sample_id, snf, family_id ->
            family_id != null && family_id != "0"
        }
        .map { sample_id, snf, family_id -> [family_id, snf] }  // Key by family_id
        .groupTuple()  // Groups by first element (family_id)
        .filter { family_id, snf_files_list -> snf_files_list.size() == 3 }
        .map { family_id, snf_files_list -> [[id: family_id], snf_files_list] }

    // Run SNIFFLES_TRIO for combined calling
    SNIFFLES_TRIO(
        ch_snf_for_trio,
        fasta
    )

    ch_trio_sv_vcf = SNIFFLES_TRIO.out.vcf
        .map { meta, vcf -> [meta + [variant_type: 'sv'], vcf] }

        RTG_MENDELIAN_SV(
            ch_trio_sv_vcf,
            ch_sdf,
            ch_ped,

        )
        ch_mendelian_vcf = RTG_MENDELIAN_SV.out.vcf
        ch_versions = ch_versions.mix(RTG_MENDELIAN_SV.out.versions)



        RTG_NONMENDELIAN_SV(
            ch_trio_sv_vcf,
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
