include { GLNEXUS                            } from '../../../modules/nf-core/glnexus/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_GLNEXUS } from '../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_GLNEXUS } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow JOINT_GENOTYPE_SNV {

    take:
    ch_gvcfs        // channel: [ val(meta), path(gvcf) ]
    ch_samplesheet  // channel: [ val(meta), val(data) ] - for family info
    ch_bed          // channel: [ val(meta), path(bed) ] or empty

    main:
    ch_versions = Channel.empty()

    // Extract family info from samplesheet
    ch_family_info = ch_samplesheet
        .map { meta, data -> [data.id, data.family_id] }
        .filter { sample_id, family_id ->
            family_id != null && family_id != "0" && family_id.trim() != ""
        }

    // Group GVCFs by family
    ch_gvcf_with_family = ch_gvcfs
        .map { meta, gvcf -> [meta.id, gvcf] }
        .join(ch_family_info, by: 0)
        .map { sample_id, gvcf, family_id -> [family_id, gvcf] }
        .groupTuple()
        .filter { family_id, gvcfs -> gvcfs.size() == 3 }  // Only trios
        .map { family_id, gvcfs ->
            [[id: family_id], gvcfs, []]  // [family_meta, gvcfs, empty custom_config]
        }

    // Run GLnexus
    GLNEXUS(
        ch_gvcf_with_family,
        ch_bed
    )

    // Index the BCF file first
    BCFTOOLS_INDEX_GLNEXUS(
        GLNEXUS.out.bcf
    )

    // Combine BCF with its index
    ch_bcf_with_index = GLNEXUS.out.bcf
        .join(BCFTOOLS_INDEX_GLNEXUS.out.csi, by: 0)

    // Convert BCF to VCF.GZ
    BCFTOOLS_VIEW_GLNEXUS(
        ch_bcf_with_index,  // Now includes [meta, bcf, index]
        [],  // regions
        [],  // targets
        []   // samples
    )

    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_GLNEXUS.out.versions)
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    emit:
    vcf = BCFTOOLS_VIEW_GLNEXUS.out.vcf
    versions = ch_versions
}
