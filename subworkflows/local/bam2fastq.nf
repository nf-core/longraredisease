// Convert merged BAM files to extract unmapped/other reads using samtools fastq with methylation tags

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'

workflow bam2fastq_subworkflow {

    take:
    ch_bam_files    // channel: [ meta, [ bam1, bam2, ... ] ]
    ch_fasta        // channel: [ meta, fasta ] (optional)
    ch_fai          // channel: [ meta, fai ] (optional)
    ch_gzi

    main:
    ch_versions = Channel.empty()

    ch_bam_files
        .branch { meta, bam_files ->
            multiple_bams: meta.is_multiple == true
                return [meta, bam_files]
            single_bam: meta.is_multiple == false
                def single_bam = bam_files instanceof List ? bam_files[0] : bam_files
                def clean_meta = [id: meta.id]
                return [clean_meta, single_bam]
        }
        .set { branched_bams }

    SAMTOOLS_MERGE (
        branched_bams.multiple_bams,
        ch_fasta,
        ch_fai,
        ch_gzi
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    ch_bams_for_fastq = SAMTOOLS_MERGE.out.bam
        .mix(branched_bams.single_bam)

    SAMTOOLS_FASTQ(
        ch_bams_for_fastq,
        false  // interleave parameter set to false
    )

    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)



emit:
    fastq       = SAMTOOLS_FASTQ.out.fastq       // channel: [meta, [fastq_1, fastq_2]] - paired-end files
    interleaved = SAMTOOLS_FASTQ.out.interleaved // channel: [meta, interleaved.fastq] - interleaved file
    singleton   = SAMTOOLS_FASTQ.out.singleton   // channel: [meta, singleton.fastq.gz] - singleton reads
    other       = SAMTOOLS_FASTQ.out.other       // channel: [meta, other.fastq.gz] - unmapped/other reads            // channel: [versions.yml]
    versions    = ch_versions                    // channel: [versions.yml] - versions of tools used in the workflow
}



