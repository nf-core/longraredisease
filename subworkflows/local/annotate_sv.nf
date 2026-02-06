include { ANNOTSV_INSTALLANNOTATIONS          } from '../../modules/nf-core/annotsv/installannotations/main.nf'
include { ANNOTSV_ANNOTSV as ANNOTSV          } from '../../modules/nf-core/annotsv/annotsv/main.nf'
include { UNTAR as UNTAR_ANNOTSV              } from '../../modules/nf-core/untar/main'

workflow annotate_sv_subworkflow {

    take:
    ch_sv_vcf       // channel: [ val(meta), path(sv_vcf) ]
    annotsv_annotations // path to pre-downloaded AnnotSV annotations tar.gz file; if not provided, annotations will be installed automatically
    candidate_genes
    false_positive_snv
    gene_transcripts

    main:
    ch_candidate_genes = candidate_genes ? Channel.fromPath(candidate_genes).map{[[id:"candidate_genes"], it]}.collect() : Channel.value([[id:"empty"], []])
    ch_false_positive_snv = false_positive_snv ? Channel.fromPath(false_positive_snv).map{[[id:"false_positive"], it]}.collect() : Channel.value([[id:"empty"], []])
    ch_gene_transcripts = gene_transcripts ? Channel.fromPath(gene_transcripts).map{[[id:"transcripts"], it]}.collect() : Channel.value([[id:"empty"], []])
    if(!annotsv_annotations) {
        ANNOTSV_INSTALLANNOTATIONS()
        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { [[id:"annotsv"], it] }
            .collect()
            .set { ch_annotsv_annotations }
    } else {
        ch_annotsv_annotations_input = Channel.fromPath(annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        if(annotsv_annotations.endsWith(".tar.gz")){
            UNTAR_ANNOTSV(ch_annotsv_annotations_input)
            UNTAR_ANNOTSV.out.untar
                .collect()
                .set { ch_annotsv_annotations }
        } else {
            ch_annotsv_annotations = Channel.fromPath(annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        }
    }

    ANNOTSV (
        ch_sv_vcf,
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    emit:
    vcf = ANNOTSV.out.vcf
    tsv = ANNOTSV.out.tsv
}
