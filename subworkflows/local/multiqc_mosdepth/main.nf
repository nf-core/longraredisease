/*
 * MultiQC subworkflow for mosdepth outputs
 * Collects mosdepth QC files and generates a MultiQC report
 */

include { MULTIQC_MOSDEPTH  } from '../../../modules/local/multiqc_mosdepth/main.nf'

workflow MULTIQC_MOSDEPTH_SUBWORKFLOW {

    take:
    mosdepth_files    // channel: [ path(file1), path(file2), ... ]


    main:
    ch_versions = channel.empty()



    // Run MultiQC
    MULTIQC_MOSDEPTH(
        mosdepth_files
    )

    ch_versions = ch_versions.mix(MULTIQC_MOSDEPTH.out.versions)

    emit:
    report   = MULTIQC_MOSDEPTH.out.report    // path: *.html
    data     = MULTIQC_MOSDEPTH.out.data      // path: *_data
    versions = ch_versions           // channel: versions
}
