// Run sniffles SV calling
include { SNIFFLES                            } from '../../modules/nf-core/sniffles/main.nf'
// Run svim SV calling
include { SVIM                                } from '../../modules/local/svim/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SVIM } from '../../modules/nf-core/bcftools/sort/main.nf'

// Run cutesv SV calling
include { CUTESV                                } from '../../modules/nf-core/cutesv/main.nf'
include { RE2SUPPORT                            } from '../../modules/local/normalize_cutesv/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CUTESV } from '../../modules/nf-core/bcftools/sort/main.nf'

workflow call_sv {

    take:
    input                    // tuple(val(meta), path(bam), path(bai))
    fasta                    // tuple(val(meta), path(fasta))
    tandem_file              // tuple(val(meta), path(bed))
    vcf_output               // val(true)
    snf_output               // val(true)


    main:
    ch_versions = Channel.empty()

    // Run all SV callers
    SNIFFLES(input, fasta, tandem_file, vcf_output, snf_output)
    SVIM(input, fasta)
    CUTESV(input, fasta)

    RE2SUPPORT(CUTESV.out.vcf)

    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(SVIM.out.versions)
    ch_versions = ch_versions.mix(CUTESV.out.versions)

    // Sort and compress VCF files from SVIM and CUTESV
    BCFTOOLS_SORT_SVIM(SVIM.out.vcf)
    BCFTOOLS_SORT_CUTESV(RE2SUPPORT.out.vcf)


    ch_versions = ch_versions.mix(BCFTOOLS_SORT_SVIM.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_CUTESV.out.versions)


    ch_sniffles_vcf_tbi = SNIFFLES.out.vcf.join(SNIFFLES.out.tbi, by: 0)
    ch_svim_vcf_tbi = BCFTOOLS_SORT_SVIM.out.vcf.join(BCFTOOLS_SORT_SVIM.out.tbi, by: 0)
    ch_cutesv_vcf_tbi = BCFTOOLS_SORT_CUTESV.out.vcf.join(BCFTOOLS_SORT_CUTESV.out.tbi, by: 0)


    emit:
    sniffles_vcf_tbi   = ch_sniffles_vcf_tbi   // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    sniffles_snf       = SNIFFLES.out.snf       // channel: [ meta, snf ]
    svim_vcf_tbi       = ch_svim_vcf_tbi       // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    cutesv_vcf_tbi     = ch_cutesv_vcf_tbi     // channel: [ meta, vcf.gz, vcf.gz.tbi ]

    sniffles_vcf     = SNIFFLES.out.vcf     // channel: [ meta, vcf.gz ]
    svim_vcf         = BCFTOOLS_SORT_SVIM.out.vcf         // channel: [ meta, vcf.gz ]
    cutesv_vcf       = BCFTOOLS_SORT_CUTESV.out.vcf       // channel: [ meta, vcf.gz ]

    // Version information
    versions             = ch_versions
}
