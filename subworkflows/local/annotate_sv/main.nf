
include { ANNOTSV_INSTALLANNOTATIONS          } from '../../../modules/nf-core/annotsv/installannotations/main.nf'
include { UNTAR as UNTAR_ANNOTSV              } from '../../../modules/nf-core/untar/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_SNIFFLES         } from '../../../modules/nf-core/annotsv/annotsv/main.nf'
include { ANNOTSV_ANNOTSV as ANNOTSV_SVIM         } from '../../../modules/nf-core/annotsv/annotsv/main.nf'

workflow ANNOTATE_SV {

    take:
    ch_samplesheet
    ch_sv_vcf
    ch_snv_vcf
    candidate_genes
    false_positive_snv
    gene_transcripts

    main:

    if(!params.annotsv_annotations) {
        ANNOTSV_INSTALLANNOTATIONS()
        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { [[id:"annotsv"], it] }
            .collect()
            .set { ch_annotsv_annotations }
    } else {
        ch_annotsv_annotations_input = Channel.fromPath(params.annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        if(params.annotsv_annotations.endsWith(".tar.gz")){
            UNTAR_ANNOTSV(ch_annotsv_annotations_input)
            UNTAR_ANNOTSV.out.untar
                .collect()
                .set { ch_annotsv_annotations }
        } else {
            ch_annotsv_annotations = Channel.fromPath(params.annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        }
    }

    ch_candidate_genes = candidate_genes ? Channel.fromPath(candidate_genes).map{[[id:"candidate_genes"], it]}.collect() : Channel.value([[id:"empty"], []])
    ch_false_positive_snv = false_positive_snv ? Channel.fromPath(false_positive_snv).map{[[id:"false_positive"], it]}.collect() : Channel.value([[id:"empty"], []])
    ch_gene_transcripts = gene_transcripts ? Channel.fromPath(gene_transcripts).map{[[id:"transcripts"], it]}.collect() : Channel.value([[id:"empty"], []])

    ch_annotate_input = ch_sv_vcf
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .join(
            ch_samplesheet.map { meta, data -> [meta.id, data.hpo_terms] },
            by: 0
        )
        .join(
            ch_snv_vcf.map { meta, vcf -> [meta.id, vcf] },
            by: 0
        )
        .map { id, meta, sv_vcf, hpo, snv_vcf ->
            def clean_hpo = (hpo && hpo.trim()) ? hpo : null
            [meta + [hpo_terms: clean_hpo], sv_vcf, [], snv_vcf]
        }

    ANNOTSV_SNIFFLES (
        ch_annotate_input,
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    if (params.run_svim){
        ANNOTSV_SVIM (ch_annotate_input,
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts)
    }

    emit:
    vcf = ANNOTSV_SNIFFLES.out.vcf
    tsv = ANNOTSV_SNIFFLES.out.tsv
    versions = ANNOTSV_SNIFFLES.out.versions
}
