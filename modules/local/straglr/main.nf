process STRAGLR {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/nourmahfel1/str-ont"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    path(bed_file)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    straglr-genotype \\
        ${bam} \\
        ${reference} \\
        --vcf ${prefix}.vcf \\
        --loci ${bed_file} \\
        --sample ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        straglr: \$(straglr-genotype --version 2>&1 | head -n1 | sed 's/.*straglr //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        straglr: \$(straglr-genotype --version 2>&1 | head -n1 | sed 's/.*straglr //')
    END_VERSIONS
    """
}
