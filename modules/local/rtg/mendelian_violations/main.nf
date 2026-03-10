process RTG_NONMENDELIAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(multisample_vcf)
    tuple val(meta2), path(sdf)
    tuple val(meta3), path(ped_file)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.versions"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rtg mendelian \\
        -i ${multisample_vcf} \\
        -t ${sdf} \\
        --pedigree ${ped_file} \\
        --output-inconsistent ${prefix}.vcf.gz

    cat <<-END_VERSIONS > ${prefix}.versions
    "${task.process}":
        rtg-tools: \$(rtg --version 2>&1 | head -1 | sed 's/.*rtg //g' | sed 's/[^0-9.].*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
