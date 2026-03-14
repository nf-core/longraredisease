include { SNPEFF_DOWNLOAD } from '../../../modules/nf-core/snpeff/download/main.nf'
include { SNPEFF_SNPEFF   } from '../../../modules/nf-core/snpeff/snpeff/main.nf'



workflow ANNOTATE_SNV {
    take:
    ch_vcf
    snpeff_db    // val: snpeff database name (e.g., 'GRCh38.105')

    main:
    ch_versions = channel.topic('versions')

    ch_db_meta = ch_vcf
        .map { meta, vcf -> [meta, snpeff_db] }.first()


    SNPEFF_DOWNLOAD(
        ch_db_meta
    )

    ch_vcf_with_cache = ch_vcf
        .combine(SNPEFF_DOWNLOAD.out.cache.map { meta, cache -> cache })
        .map { meta, vcf, cache -> [meta, vcf, [id: 'cache'], cache] }

    // Annotate VCF
    SNPEFF_SNPEFF(
        ch_vcf_with_cache.map { meta, vcf, meta2, cache -> [meta, vcf] },
        snpeff_db,
        ch_vcf_with_cache.map { meta, vcf, meta2, cache -> [meta2, cache] }
    )


    emit:
    vcf          = SNPEFF_SNPEFF.out.vcf           // channel: [ val(meta), path(vcf) ]
    report       = SNPEFF_SNPEFF.out.report        // channel: [ val(meta), path(csv) ]
    summary_html = SNPEFF_SNPEFF.out.summary_html  // channel: [ val(meta), path(html) ]
    genes_txt    = SNPEFF_SNPEFF.out.genes_txt     // channel: [ val(meta), path(txt) ]
    cache        = SNPEFF_DOWNLOAD.out.cache       // channel: [ val(meta), path(cache) ]
    versions     = ch_versions                      // channel: versions
}
