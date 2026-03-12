process MULTIQC_MOSDEPTH {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b'

    input:
    tuple val(meta), path(input_files, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}_multiqc.html"), emit: report
    tuple val(meta), path("${meta.id}_multiqc_data"), emit: data, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc \\
        -v \\
        -o . \\
        -n ${meta.id}_multiqc.html \\
        ${args} \\
        input/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/.* //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.id}_multiqc_data
    touch ${meta.id}_multiqc.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/.* //g')
    END_VERSIONS
    """
}
