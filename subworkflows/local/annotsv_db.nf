include { ANNOTSV_INSTALLANNOTATIONS          } from '../../modules/nf-core/annotsv/installannotations/main.nf'
include { UNTAR as UNTAR_ANNOTSV              } from '../../modules/nf-core/untar/main'

workflow annotsv_db {  // Fixed typo: annptsv_dv -> annotsv_db

    take:
    annotsv_annotations // path to pre-downloaded AnnotSV annotations tar.gz file

    main:

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

    emit:
    db = ch_annotsv_annotations

}
