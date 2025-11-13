include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_SV } from '../../modules/nf-core/bcftools/view/main.nf'
include { DOWNSAMPLE_SV } from '../../modules/local/downsample_sv/main.nf'

workflow filter_sv {

take:

    ch_vcf_tbi        // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    target_bed        // channel: [ meta, bed ] (optional)
    downsample_sv    // val: true/false
    mosdepth_summary  // channel: [ meta, mosdepth_summary.txt ]
    mosdepth_bed      // channel: [ meta, mosdepth_quantized.bed ] (optional)
    chromosome_codes  // val: list of chromosome codes
    min_read_support  // val: minimum read support
    min_read_support_limit // val: minimum read support limit


main:

ch_versions = Channel.empty()
target_bed_ch = target_bed ?: Channel.value([])


BCFTOOLS_FILTER_SV (
    ch_vcf_tbi,
    Channel.value([]),
    target_bed_ch,
    Channel.value([])
    )

ch_sv = BCFTOOLS_FILTER_SV.out.vcf.join(BCFTOOLS_FILTER_SV.out.tbi, by: 0)
ch_versions = ch_versions.mix(BCFTOOLS_FILTER_SV.out.versions)


    if (downsample_sv) {
        // Clean metadata for join (remove caller info temporarily)
        ch_sv_clean_meta = ch_sv.map { meta, vcf, tbi ->
            def clean_meta = [id: meta.id]  // Keep only sample ID
            def original_meta = meta        // Store original metadata
            [clean_meta, vcf, tbi, original_meta]
        }

        // Join with mosdepth data using clean metadata
        ch_sv_mosdepth_bed = ch_sv_clean_meta
            .map { clean_meta, vcf, tbi, original_meta -> [clean_meta, vcf, tbi, original_meta] }
            .join(mosdepth_summary, by: 0)  // Join by clean sample ID
            .join(mosdepth_bed, by: 0)      // Join by clean sample ID
            .map { clean_meta, vcf, tbi, original_meta, summary, bed ->
                [original_meta, vcf, tbi, summary, bed]  // Restore original metadata with caller info
            }

        DOWNSAMPLE_SV(
            ch_sv_mosdepth_bed.map { meta, vcf, tbi, summary, bed -> tuple(meta, vcf, tbi) },
            ch_sv_mosdepth_bed.map { meta, vcf, tbi, summary, bed -> tuple(meta, summary) },
            ch_sv_mosdepth_bed.map { meta, vcf, tbi, summary, bed -> tuple(meta, bed) },
            chromosome_codes,
            min_read_support,
            min_read_support_limit
        )
        ch_sv = DOWNSAMPLE_SV.out.filterbycov_vcf
        ch_versions = ch_versions.mix(DOWNSAMPLE_SV.out.versions)
    }

emit:
ch_vcf_tbi = ch_sv   // channel: [ meta, vcf.gz, vcf.gz.tbi ]
versions = ch_versions
}
