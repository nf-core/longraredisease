// This workflow is for clair3

include { CLAIR3 } from '../../modules/nf-core/clair3/main.nf'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../modules/nf-core/deepvariant/rundeepvariant/main.nf'
include { DEEPVARIANT_VCFSTATSREPORT } from '../../modules/nf-core/deepvariant/vcfstatsreport/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_CLAIR3 } from '../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_DEEPVARIANT } from '../../modules/nf-core/bcftools/view/main.nf'

workflow call_snv {
    take:
    ch_input_bam            // channel: tuple(val(meta), path(bam), path(bai))
    fasta                   // channel: tuple(val(meta2), path(fasta))
    fai                     // channel: tuple(val(meta3), path(fai)) - optional
    run_deepvariant         // boolean
    ch_input_deepvariant

    main:
    ch_versions = Channel.empty()

    ch_vcf = Channel.empty()
    ch_tbi = Channel.empty()
    ch_gvcf = Channel.empty()
    ch_gtbi = Channel.empty()
    
    // Model and platform presets
    model_presets = [
        ont: 'ont',
        pacbio: 'hifi',
        hifi: 'hifi'
    ]
    
    platform_presets = [
        ont: 'ont',
        pacbio: 'hifi',
        hifi: 'hifi'
    ]
    
    // Resolve model and platform
    resolved_model = params.clair3_model ?: model_presets[params.sequencing_platform] ?: 'ont'
    resolved_platform = platform_presets[params.sequencing_platform] ?: 'ont'
    
    // Prepare input channel with resolved values
    ch_input_clair3 = ch_input_bam.map { meta, bam, bai ->
        tuple(
            meta,
            bam,
            bai,
            resolved_model,      // packaged_model
            [],                  // user_model (empty for packaged models)
            resolved_platform    // platform
        )
    }
    
    CLAIR3(
        ch_input_clair3,    // tuple(meta, bam, bai, packaged_model, user_model, platform)
        fasta,              // tuple(meta2, fasta)
        fai                 // tuple(meta3, fai)
    )
    
    // Handle CLAIR3 filtering
    if (params.filter_pass_snv) {
        ch_clair3_vcf = CLAIR3.out.vcf
            .join(CLAIR3.out.tbi, by: 0)

        BCFTOOLS_FILTER_CLAIR3(
            ch_clair3_vcf,      // tuple(meta, vcf, tbi)
            Channel.value([]),   // empty channel for samples
            Channel.value([]),   // empty channel for regions
            Channel.value([])    // empty channel for filters
        )

        ch_vcf = BCFTOOLS_FILTER_CLAIR3.out.vcf
        ch_tbi = BCFTOOLS_FILTER_CLAIR3.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_CLAIR3.out.versions)

    } else {
        ch_vcf = CLAIR3.out.vcf
        ch_tbi = CLAIR3.out.tbi
    }

    ch_gvcf = CLAIR3.out.gvcf
    ch_gtbi = CLAIR3.out.gtbi
    ch_versions = ch_versions.mix(CLAIR3.out.versions)


    if (run_deepvariant) {
        DEEPVARIANT_RUNDEEPVARIANT(
            ch_input_deepvariant,
            fasta,
            fai,
            [[:], []],
            [[:], []]
        )

        if (params.filter_pass_snv) {
            ch_deepvariant_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
                .join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index, by: 0)

            BCFTOOLS_FILTER_DEEPVARIANT(
                ch_deepvariant_vcf,
                Channel.value([]),
                Channel.value([]),
                Channel.value([])
            )

            ch_vcf_deepvariant = BCFTOOLS_FILTER_DEEPVARIANT.out.vcf
            ch_tbi_deepvariant = BCFTOOLS_FILTER_DEEPVARIANT.out.tbi
            ch_versions = ch_versions.mix(BCFTOOLS_FILTER_DEEPVARIANT.out.versions)

        } else {
            ch_vcf_deepvariant = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
            ch_tbi_deepvariant = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index
        }
        
        ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)

        if (params.deepvariant_runtime_report){
            DEEPVARIANT_VCFSTATSREPORT(DEEPVARIANT_RUNDEEPVARIANT.out.vcf)
            html_report = DEEPVARIANT_VCFSTATSREPORT.out.report
        }
        else {
            html_report = Channel.empty
        }
        
        } 
        
        else {
        ch_vcf_deepvariant = Channel.empty()
        ch_tbi_deepvariant = Channel.empty()
        
        }

    emit:
    vcf             = ch_vcf      
    tbi             = ch_tbi      
    gvcf            = ch_gvcf        
    gtbi            = ch_gtbi 
    phased_vcf      = CLAIR3.output.phased_vcf
    phased_tbi      = CLAIR3.output.phased_tbi
    deepvariant_vcf = ch_vcf_deepvariant
    deepvariant_tbi = ch_tbi_deepvariant
    deepvariant_report  = html_report       
    versions        = ch_versions
}