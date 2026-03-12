
include { ANNOTSV_ANNOTSV as ANNOTSV          } from '../../modules/nf-core/annotsv/annotsv/main.nf'


workflow annotate_sv {

    take:
    ch_sv_vcf           // channel: [ val(meta), path(sv_vcf) ] - NO INDEX!
    ch_hpo_terms        // channel: [ val(meta), hpo_terms ]
    ch_snv_vcf          // channel: [ val(meta), path(snv_vcf) ]
    ch_annotsv_annotations // path to pre-downloaded AnnotSV annotations tar.gz file
    candidate_genes
    false_positive_snv
    gene_transcripts

    main:
    ch_candidate_genes = candidate_genes ? Channel.fromPath(candidate_genes).map{[[id:"candidate_genes"], it]}.collect() : Channel.value([[id:"empty"], []])
    ch_false_positive_snv = false_positive_snv ? Channel.fromPath(false_positive_snv).map{[[id:"false_positive"], it]}.collect() : Channel.value([[id:"empty"], []])
    ch_gene_transcripts = gene_transcripts ? Channel.fromPath(gene_transcripts).map{[[id:"transcripts"], it]}.collect() : Channel.value([[id:"empty"], []])

    // Join SV VCF with HPO terms by meta.id
    // Clean HPO terms - set to null if empty string
    ch_sv_with_hpo = ch_sv_vcf
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .join(
            ch_hpo_terms.map { meta, hpo -> [meta.id, hpo] },
            by: 0
        )
        .map { id, meta, sv_vcf, hpo ->
            def clean_hpo = (hpo && hpo.trim()) ? hpo : null
            [meta + [hpo_terms: clean_hpo], sv_vcf]
        }

    // Join with SNV VCF to get candidate small variants
    ch_sv_with_snv = ch_sv_with_hpo
        .map { meta, sv_vcf -> [meta.id, meta, sv_vcf] }
        .join(
            ch_snv_vcf.map { meta, vcf -> [meta.id, vcf] },
            by: 0
        )
        .map { id, meta, sv_vcf, snv_vcf ->
            [meta, sv_vcf, [], snv_vcf]  // Empty list for sv_vcf_index
        }

    ANNOTSV (
        ch_sv_with_snv,
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    emit:
    vcf = ANNOTSV.out.vcf
    tsv = ANNOTSV.out.tsv
    versions = ANNOTSV.out.versions
}
