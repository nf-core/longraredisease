include {
    SAMTOOLS_MERGE
} from '../../modules/nf-core/samtools/merge/main'
include {
    SAMTOOLS_FASTQ
} from '../../modules/nf-core/samtools/fastq/main'
include {
    SAMTOOLS_VIEW
} from '../../modules/nf-core/samtools/view/main'

workflow bam2fastq_subworkflow {
    take:
    ch_bam_files    // channel: [ meta, [ bam1, bam2, ... ] ]
    ch_fasta        // channel: [ meta, fasta ] (optional)
    ch_fai          // channel: [ meta, fai ] (optional)
    ch_gzi          // channel: [ meta, gzi ] (optional)

    main:
    ch_versions = Channel.empty()

    // Run samtools merge
    SAMTOOLS_MERGE (
        ch_bam_files,
        ch_fasta,
        ch_fai,
        ch_gzi
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // Prepare input for SAMTOOLS_VIEW - add empty index to merged BAM
    ch_merged_bam_with_index = SAMTOOLS_MERGE.out.bam.map { meta, bam ->
        [meta, bam, []]  // Add empty index as third element
    }

    // Extract unmapped reads and convert to CRAM
    SAMTOOLS_VIEW (
        ch_merged_bam_with_index,         // tuple val(meta), path(input), path(index)
        ch_fasta,                         // tuple val(meta2), path(fasta)
        [[:], []],                  // path qname (no qname file)
        "crai"                           // val index_format
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    // Convert merged BAM to FASTQ
    SAMTOOLS_FASTQ(
        SAMTOOLS_MERGE.out.bam,
        false  // interleave parameter set to false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    emit:
    // Original outputs
    fastq = SAMTOOLS_FASTQ.out.fastq           // channel: [meta, [fastq_1, fastq_2]]
    interleaved = SAMTOOLS_FASTQ.out.interleaved // channel: [meta, interleaved.fastq]
    singleton = SAMTOOLS_FASTQ.out.singleton   // channel: [meta, singleton.fastq.gz]
    other = SAMTOOLS_FASTQ.out.other           // channel: [meta, other.fastq.gz]

    // New unmapped CRAM output
    unmapped_cram = SAMTOOLS_VIEW.out.cram     // channel: [meta, unmapped.cram]
    unmapped_crai = SAMTOOLS_VIEW.out.crai     // channel: [meta, unmapped.cram.crai]

    versions = ch_versions                      // channel: [versions.yml]
}
