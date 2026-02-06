// Run sniffles SV calling
include { SNIFFLES                                   } from '../../modules/nf-core/sniffles/main.nf'
include { GUNZIP as GUNZIP_SNIFFLES_PLOT             } from '../../modules/nf-core/gunzip/main.nf'
include { SV_PLOT                                    } from '../../modules/local/generate_sv_plots/main.nf'
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
    single_caller            // val(boolean) - whether to run single caller mode
    sv_caller                // val(string) - which caller to use in single mode
    generate_sniffles_plots           // val(boolean) - whether to generate SV plots

    main:
    ch_versions = Channel.empty()

    // Initialize empty channels for all callers
    ch_sniffles_vcf = Channel.empty()
    ch_sniffles_tbi = Channel.empty()
    ch_sniffles_snf = Channel.empty()
    ch_svim_vcf = Channel.empty()
    ch_svim_tbi = Channel.empty()
    ch_cutesv_vcf = Channel.empty()
    ch_cutesv_tbi = Channel.empty()
    ch_sniffles_plots = Channel.empty()

    if (single_caller) {
        // Single caller mode
        if (sv_caller == 'sniffles') {
            SNIFFLES(input, fasta, tandem_file, vcf_output, snf_output)
            ch_sniffles_vcf = SNIFFLES.out.vcf
            ch_sniffles_tbi = SNIFFLES.out.tbi
            ch_sniffles_snf = SNIFFLES.out.snf
            ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        } else if (sv_caller == 'svim') {
            SVIM(input, fasta)
            BCFTOOLS_SORT_SVIM(SVIM.out.vcf)
            ch_svim_vcf = BCFTOOLS_SORT_SVIM.out.vcf
            ch_svim_tbi = BCFTOOLS_SORT_SVIM.out.tbi
            ch_versions = ch_versions.mix(SVIM.out.versions)
            ch_versions = ch_versions.mix(BCFTOOLS_SORT_SVIM.out.versions)

        } else if (sv_caller == 'cutesv') {
            CUTESV(input, fasta)
            RE2SUPPORT(CUTESV.out.vcf)
            BCFTOOLS_SORT_CUTESV(RE2SUPPORT.out.vcf)
            ch_cutesv_vcf = BCFTOOLS_SORT_CUTESV.out.vcf
            ch_cutesv_tbi = BCFTOOLS_SORT_CUTESV.out.tbi
            ch_versions = ch_versions.mix(CUTESV.out.versions)
            ch_versions = ch_versions.mix(BCFTOOLS_SORT_CUTESV.out.versions)

        } else {
            error "Invalid sv_caller specified: ${sv_caller}. Valid options are: sniffles, svim, cutesv"
        }

    } else {
        // Multi-caller mode (original behavior)
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

        ch_sniffles_vcf = SNIFFLES.out.vcf
        ch_sniffles_tbi = SNIFFLES.out.tbi
        ch_sniffles_snf = SNIFFLES.out.snf
        ch_svim_vcf = BCFTOOLS_SORT_SVIM.out.vcf
        ch_svim_tbi = BCFTOOLS_SORT_SVIM.out.tbi
        ch_cutesv_vcf = BCFTOOLS_SORT_CUTESV.out.vcf
        ch_cutesv_tbi = BCFTOOLS_SORT_CUTESV.out.tbi
    }

    if (generate_sniffles_plots && (!single_caller || sv_caller == 'sniffles')) {
        GUNZIP_SNIFFLES_PLOT(ch_sniffles_vcf)
        SV_PLOT(GUNZIP_SNIFFLES_PLOT.out.gunzip)
        ch_sniffles_plots = SV_PLOT.out.plot_dir
        ch_versions = ch_versions.mix(SV_PLOT.out.versions)
    }

    // Combine VCF and TBI channels
    ch_sniffles_vcf_tbi = ch_sniffles_vcf.join(ch_sniffles_tbi, by: 0, remainder: true)
        .filter { meta, vcf, tbi -> vcf != null }
    ch_svim_vcf_tbi = ch_svim_vcf.join(ch_svim_tbi, by: 0, remainder: true)
        .filter { meta, vcf, tbi -> vcf != null }
    ch_cutesv_vcf_tbi = ch_cutesv_vcf.join(ch_cutesv_tbi, by: 0, remainder: true)
        .filter { meta, vcf, tbi -> vcf != null }

    emit:
    sniffles_vcf_tbi   = ch_sniffles_vcf_tbi   // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    sniffles_snf       = ch_sniffles_snf       // channel: [ meta, snf ]
    svim_vcf_tbi       = ch_svim_vcf_tbi       // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    cutesv_vcf_tbi     = ch_cutesv_vcf_tbi     // channel: [ meta, vcf.gz, vcf.gz.tbi ]

    sniffles_vcf     = ch_sniffles_vcf     // channel: [ meta, vcf.gz ]
    svim_vcf         = ch_svim_vcf         // channel: [ meta, vcf.gz ]
    cutesv_vcf       = ch_cutesv_vcf       // channel: [ meta, vcf.gz ]

    // Version information
    versions             = ch_versions
}
