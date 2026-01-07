include { STRANGER   } from '../../modules/nf-core/stranger/main'


workflow annotate_str {
    take:
    vcf      // channel: [ val(sample_id), path(bam) ]
    variant_catalogue  // channel: path(variant_catalogue)

    main:
    STRANGER(
        vcf,
        variant_catalogue)

    emit:
    vcf      = STRANGER.out.vcf      // channel: [ val(meta), path(vcf) ]
    tbi      = STRANGER.out.tbi      // channel: [ val
    versions = STRANGER.out.versions_stranger // channel: path(versions.yml)

    }
