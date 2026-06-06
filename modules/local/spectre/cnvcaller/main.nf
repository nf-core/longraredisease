process SPECTRE_CNVCALLER {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/ont-spectre:0.3.2--adfae189059be3d9"

    input:
    tuple val(meta), path(mosdepth_summary), path(mosdepth_regions_bed), path(mosdepth_regions_csi), path(vcf)
    tuple val(meta2), path(fasta)
    path(metadata_file)
    path(blacklist)
    val bin_size

    output:
    tuple val(meta), path("*.vcf")                       , emit: vcf
    tuple val(meta), path("*.bed")                       , emit: bed
    tuple val(meta), path("*.spc.gz")                    , emit: spc
    tuple val(meta), path("predicted_karyotype.txt")     , emit: txt
    tuple val(meta), path("windows_stats", type: 'dir')  , emit: winstats
    tuple val("${task.process}"), val('spectre'), eval('spectre version 2>&1 | grep -oP "Spectre version: \\K[0-9.]+"'), emit: versions_spectre, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def coverage_dir = "mosdepth_${prefix}"

    """
    mkdir -p ${coverage_dir}

    ln -s \$(realpath ${mosdepth_summary}) ${coverage_dir}/${meta.id}.mosdepth.summary.txt
    ln -s \$(realpath ${mosdepth_regions_bed}) ${coverage_dir}/${meta.id}.regions.bed.gz
    ln -s \$(realpath ${mosdepth_regions_csi}) ${coverage_dir}/${meta.id}.regions.bed.gz.csi

    # Now run spectre with the directory - use bin_size input parameter, not params
    spectre CNVCaller \\
        --coverage ${coverage_dir} \\
        --snv ${vcf} \\
        --reference ${fasta} \\
        --output-dir . \\
        --metadata ${metadata_file} \\
        --blacklist ${blacklist} \\
        --bin-size ${bin_size} \\
        ${args}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.bed
    touch ${prefix}.spc.gz
    touch predicted_karyotype.txt
    mkdir windows_stats

    """
}
