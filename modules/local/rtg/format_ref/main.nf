process RTG_FORMAT_REF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.sdf"), emit: sdf
    tuple val(meta), path("*.versions"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }

    """
    rtg RTG_MEM=${avail_mem} format \\
        ${args} \\
        -o ${prefix}.sdf \\
        ${fasta}

    cat <<-END_VERSIONS > ${prefix}.versions
    "${task.process}":
        rtg-tools: \$(rtg --version 2>&1 | head -1 | sed 's/.*rtg //g' | sed 's/[^0-9.].*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sdf
    touch versions.yml
    """
}
