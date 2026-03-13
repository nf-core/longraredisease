// Run sniffles SV calling
include { SNIFFLES                                   } from '../../modules/nf-core/sniffles/main.nf'
include { GUNZIP as GUNZIP_SNIFFLES_PLOT             } from '../../modules/nf-core/gunzip/main.nf'
include { SNIFFLES_GENERATE_PLOTS                    } from '../../modules/local/sniffles/generate_plots/main.nf'
// Run svim SV calling
include { SVIM_ALIGNMENT                      } from '../../modules/local/svim/alignment/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SVIM } from '../../modules/nf-core/bcftools/sort/main.nf'
// Run cutesv SV calling
include { CUTESV                                } from '../../modules/nf-core/cutesv/main.nf'
include { RE2SUPPORT                            } from '../../modules/local/fix_header_sv/cutesv/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CUTESV } from '../../modules/nf-core/bcftools/sort/main.nf'

workflow call_sv {

    take:
    input                    // tuple(val(meta), path(bam), path(bai))
    fasta                    // tuple(val(meta), path(fasta))
    tandem_file              // tuple(val(meta), path(bed))
    vcf_output               // val(true)
    snf_output               // val(true)
    merge_sv                 // val(boolean) - whether to prepare for merging (i.e. run all callers regardless of run_svim/run_cutesv)
    run_svim

    main:
    ch_versions = channel.empty()

    // Initialize empty channels for conditional callers
    ch_svim_vcf = channel.empty()
    ch_svim_tbi = channel.empty()
    ch_cutesv_vcf = channel.empty()
    ch_cutesv_tbi = channel.empty()

    // ========================================
    // SNIFFLES - ALWAYS RUNS
    // ========================================
    SNIFFLES(input, fasta, tandem_file, vcf_output, snf_output)
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)

    // ========================================
    // SNIFFLES PLOTS
    // ========================================
    GUNZIP_SNIFFLES_PLOT(SNIFFLES.out.vcf)
    SNIFFLES_GENERATE_PLOTS(GUNZIP_SNIFFLES_PLOT.out.gunzip)
    ch_sniffles_plots = SNIFFLES_GENERATE_PLOTS.out.plot_dir


    if (merge_sv || run_svim) {

    // ========================================
    // SVIM - CONDITIONAL
    // ========================================

        SVIM_ALIGNMENT(input, fasta)
        BCFTOOLS_SORT_SVIM(SVIM_ALIGNMENT.out.vcf)

        ch_svim_vcf = BCFTOOLS_SORT_SVIM.out.vcf
        ch_svim_tbi = BCFTOOLS_SORT_SVIM.out.tbi
        }


    // ========================================
    // CUTESV
    // ========================================
        if (merge_sv) {
            CUTESV(input, fasta)
            RE2SUPPORT(CUTESV.out.vcf)
            BCFTOOLS_SORT_CUTESV(RE2SUPPORT.out.vcf)

            ch_cutesv_vcf = BCFTOOLS_SORT_CUTESV.out.vcf
            ch_cutesv_tbi = BCFTOOLS_SORT_CUTESV.out.tbi
            ch_versions = ch_versions.mix(CUTESV.out.versions)

            }

    // ========================================
    // COMBINE VCF AND TBI CHANNELS
    // ========================================

    ch_sniffles_vcf_tbi = SNIFFLES.out.vcf.join(SNIFFLES.out.tbi, by: 0)

    ch_svim_vcf_tbi = ch_svim_vcf.join(ch_svim_tbi, by: 0, remainder: true)
        .filter { meta, vcf, tbi -> vcf != null }

    ch_cutesv_vcf_tbi = ch_cutesv_vcf.join(ch_cutesv_tbi, by: 0, remainder: true)
        .filter { meta, vcf, tbi -> vcf != null }


    emit:
    sniffles_vcf_tbi = ch_sniffles_vcf_tbi   // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    sniffles_snf     = SNIFFLES.out.snf      // channel: [ meta, snf ]
    sniffles_vcf     = SNIFFLES.out.vcf      // channel: [ meta, vcf.gz ]
    sniffles_unzipped_vcf = GUNZIP_SNIFFLES_PLOT.out.gunzip
    svim_vcf_tbi     = ch_svim_vcf_tbi       // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    svim_vcf         = ch_svim_vcf           // channel: [ meta, vcf.gz ]
    cutesv_vcf_tbi   = ch_cutesv_vcf_tbi     // channel: [ meta, vcf.gz, vcf.gz.tbi ]
    cutesv_vcf       = ch_cutesv_vcf         // channel: [ meta, vcf.gz ]
    sniffles_plots   = ch_sniffles_plots     // channel: [ meta, plot_dir ]
    versions         = ch_versions           // channel: [ versions ]
}
