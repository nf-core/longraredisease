include { JASMINESV } from '../../modules/nf-core/jasminesv/main'
include { NORMALIZE_JASMINE } from '../../modules/local/normalize_merged_sv/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_JASMINE } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_JASMINE } from '../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_JASMINE } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX_JASMINE } from '../../modules/nf-core/tabix/tabix/main'

workflow merge_sv {
    take:
    ch_vcfs         // channel: [ meta, [vcf1, vcf2, vcf3] ] - VCFs from 3 aligners
    ch_fasta        // channel: [ meta, fasta ]
    ch_fasta_fai    // channel: [ meta, fasta_fai ]
    ch_chr_norm     // channel: path(chr_norm) (optional)

    main:
    ch_versions = Channel.empty()

    // Step 1: Run JASMINESV to merge VCFs
    JASMINESV(
        ch_vcfs,
        ch_fasta,
        ch_fasta_fai,
        ch_chr_norm ?: Channel.value([])
    )
    ch_versions = ch_versions.mix(JASMINESV.out.versions)

    // Step 2: Normalize the merged VCF and create sample file
    NORMALIZE_JASMINE(JASMINESV.out.vcf)
    // ch_versions = ch_versions.mix(NORMALIZE_JASMINE.out.versions)

    // Step 3: Sort the normalized VCF
    jasmine_vcf = NORMALIZE_JASMINE.out.vcf_with_samples.map { meta, vcf, samples -> [meta, vcf] }
    BCFTOOLS_SORT_JASMINE(jasmine_vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_JASMINE.out.versions)

    // Step 4: View to ensure proper compression and indexing
    BCFTOOLS_VIEW_JASMINE(
        BCFTOOLS_SORT_JASMINE.out.vcf.join(BCFTOOLS_SORT_JASMINE.out.tbi, by: 0),
        Channel.value([]),  // regions - use Channel.value([]) instead of []
        Channel.value([]),  // targets
        Channel.value([])   // samples
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_JASMINE.out.versions)

    // Step 5: Reheader with proper sample name
    BCFTOOLS_REHEADER_JASMINE(
        BCFTOOLS_VIEW_JASMINE.out.vcf
            .join(
                NORMALIZE_JASMINE.out.vcf_with_samples.map { meta, vcf, samples -> [meta, samples] },
                by: 0
            )
            .map { meta, vcf, samples ->
                [meta, vcf, [], samples]  // [meta, vcf, header, samples]
            },
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER_JASMINE.out.versions)

    // Step 6: Index final VCF
    TABIX_JASMINE(
        BCFTOOLS_REHEADER_JASMINE.out.vcf.map { meta, vcf -> [meta, vcf] }
    )
    ch_versions = ch_versions.mix(TABIX_JASMINE.out.versions)
    emit:
    vcf      = BCFTOOLS_REHEADER_JASMINE.out.vcf  // Changed to final output
    tbi      = TABIX_JASMINE.out.tbi
    versions = ch_versions
}
