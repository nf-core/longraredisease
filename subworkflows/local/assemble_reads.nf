include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main'
include { HIFIASM   } from '../../modules/nf-core/hifiasm/main'
include { YAK_COUNT } from '../../modules/nf-core/yak/count/main'
include { GFASTATS  } from '../../modules/nf-core/gfastats/main'

// This subworkflow performs genome assembly and generates haplotype outputs from sequencing reads (organized by sample), utilizing hifiasm and gfastats tools.
// The workflow assumes each sample may contain multiple files, but each sample is associated with at most one family unit.
workflow GENOME_ASSEMBLY {

    take:
    input_reads        // channel: [ val(meta), fastqs ]
    enable_trio_analysis   //    bool: Should we use trio binning mode where possible?

    main:
    software_versions = Channel.empty()

    if (enable_trio_analysis) {
        // Initially, we categorize samples according to their familial relationships
        input_reads
            .branch { sample_meta, _fastq_files ->
                def parent_status = sample_meta.relationship in ['father', 'mother']
                parent_pairs                    : parent_status && sample_meta.has_other_parent
                offspring_with_complete_parents : sample_meta.relationship == 'child' && sample_meta.two_parents
                remaining_samples               : true
            }
            .set { categorized_samples }

        // Subsequently, files from parents of offspring with complete parent pairs require concatenation prior to yak processing
        // when multiple files exist for the same parent.
        categorized_samples.parent_pairs
            .branch { _sample_meta, fastq_list ->
                requires_concat: fastq_list.size() > 1
                single_file: fastq_list.size() == 1
            }
            .set { parent_files_for_yak }

        CAT_FASTQ (
            parent_files_for_yak.requires_concat
        )
        software_versions = software_versions.mix(CAT_FASTQ.out.versions)

        YAK_COUNT (
            CAT_FASTQ.out.reads.concat(parent_files_for_yak.single_file)
        )
        software_versions = software_versions.mix(YAK_COUNT.out.versions)

        YAK_COUNT.out.yak
            // Since a parent may have multiple offspring, and sample_meta.children contains all offspring,
            // we must generate one tuple per offspring.
            .flatMap { sample_meta, yak_file ->
                (sample_meta.children ?: []).collect { offspring_id ->
                    [offspring_id, sample_meta, yak_file]
                }
            }
            .branch { offspring_id, sample_meta, yak_file ->
                father_line: sample_meta.relationship == 'father'
                    return [ offspring_id, yak_file ]
                mother_line: sample_meta.relationship == 'mother'
                    return [ offspring_id, yak_file ]
            }
            .set { yak_results }

        // Constructs input for trio-binned assemblies (offspring with complete parent pairs)
        categorized_samples.offspring_with_complete_parents
            .map { sample_meta, fastq_files -> [ sample_meta.id, sample_meta, fastq_files ] }
            .join(yak_results.father_line)
            .join(yak_results.mother_line)
            .map { _sample_id, sample_meta, fastq_files, paternal_yak, maternal_yak ->
                [ sample_meta, fastq_files, paternal_yak, maternal_yak ]
            }
            .set { complete_parent_data }

        // Construct hifiasm input by merging non-trio binned samples with trio-binned samples.
        categorized_samples.remaining_samples
            .concat(categorized_samples.parent_pairs)
            .map { sample_meta, fastq_list ->
                [ sample_meta, fastq_list, [], [] ]
            }
            .concat(complete_parent_data)
            .multiMap { sample_meta, fastq_files, paternal_yak, maternal_yak ->
                sequence_data : [ sample_meta, fastq_files  , []          ]
                yak_data      : [ sample_meta, paternal_yak , maternal_yak ]
            }
            .set { hifiasm_input_channels }
    } else {
        input_reads
            .multiMap { sample_meta, fastq_files ->
                sequence_data : [ sample_meta, fastq_files, [] ]
                yak_data      : [ [], [], [] ]
            }
            .set { hifiasm_input_channels }
    }

    HIFIASM (
        hifiasm_input_channels.sequence_data,
        hifiasm_input_channels.yak_data,
        [[],[],[]],
        [[],[]]
    )
    software_versions = software_versions.mix(HIFIASM.out.versions)

    HIFIASM.out.hap1_contigs
        .map { sample_meta, assembly_fasta -> [ sample_meta + [ 'haplotype': 1 ], assembly_fasta ] }
        .set { paternal_gfastats_input }

    HIFIASM.out.hap2_contigs
        .map { sample_meta, assembly_fasta -> [ sample_meta + [ 'haplotype': 2 ], assembly_fasta ] }
        .set { maternal_gfastats_input }

    GFASTATS(
        paternal_gfastats_input.mix(maternal_gfastats_input),
        'fasta',
        '',
        '',
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    software_versions = software_versions.mix(GFASTATS.out.versions)

    emit:
    assembled_haplotypes = GFASTATS.out.assembly // channel: [ val(meta), path(fasta) ]
    versions = software_versions                 // channel: [ versions.yml ]
}
