// Unify multiple structural variant VCF files (SV, CNV, repeat)

include { UNIFY_GENEYX } from '../../modules/local/unifyvcf_geneyx/main'
include {GUNZIP as GUNZIP_UNIFY} from '../../modules/nf-core/gunzip/main'

workflow unify_geneyx {

    take:
    ch_sv_vcfs          // channel: [meta, [sv1.vcf, sv2.vcf, ...]] - Multiple SV VCF files
    ch_cnv_vcf          // channel: [meta, cnv.vcf] - Single CNV VCF file (optional)
    ch_repeat_vcf       // channel: [meta, repeat.vcf] - Single repeat VCF file (optional)
    modify_repeats      // Boolean: whether to modify repeat calls (true for STRaglr)

    main:

    ch_versions = Channel.empty()

    // Unify all VCF files
    UNIFY_GENEYX(
        ch_sv_vcfs,
        ch_cnv_vcf,
        ch_repeat_vcf,
        modify_repeats
    )
// Just for scientists to view it
    GUNZIP_UNIFY(
        UNIFY_GENEYX.out.unified_vcf
    )

    ch_versions = ch_versions.mix(UNIFY_GENEYX.out.versions)

    emit:
    unified_vcf = GUNZIP_UNIFY.out.gunzip // channel: [meta, unified.vcf.gz]
    unified_vcf_gz = UNIFY_GENEYX.out.unified_vcf    // channel: [meta, unified.vcf]
    unified_tbi = UNIFY_GENEYX.out.unified_tbi    // channel: [meta, unified.vcf.tbi]
    versions    = ch_versions                   // channel: versions.yml
}
