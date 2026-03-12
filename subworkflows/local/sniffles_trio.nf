include { SNIFFLES_TRIO                      } from '../../modules/local/sniffles/trio/main.nf'


workflow sniffles_trio {
    take:
    snf_files     // channel: [ val(meta), path(snf) ]
    samplesheet   // channel: [ val(meta), val(data) ] - with family_id info
    fasta         // channel: path(fasta)

    main:

    // Group SNF files by family for trio calling
    ch_snf_for_trio = snf_files
        .map { meta, snf -> [meta.id, snf] }
        .join(
            samplesheet.map { meta, data -> [data.id, data.family_id] },
            by: 0
        )
        .filter { sample_id, snf, family_id ->
            family_id != null && family_id != "0"
        }
        .map { sample_id, snf, family_id -> [family_id, snf] }  // Key by family_id
        .groupTuple()  // Groups by first element (family_id)
        .filter { family_id, snf_files_list -> snf_files_list.size() == 3 }
        .map { family_id, snf_files_list -> [[id: family_id], snf_files_list] }

    // Run SNIFFLES_TRIO for combined calling
    SNIFFLES_TRIO(
        ch_snf_for_trio,
        fasta
    )



    emit:
    vcf      = SNIFFLES_TRIO.out.vcf      // channel: [ val(meta), path(vcf) ]



}
