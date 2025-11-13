process ANNOTATE_STR {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::python=3.9 bioconda::htslib=1.21"
    container "community.wave.seqera.io/library/htslib_procs_python:4c242671a1021f9c"

    input:
    tuple val(meta), path(vcf)
    path variant_catalogue

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_str.py \\
        --vcf ${vcf} \\
        --catalogue ${variant_catalogue} \\
        --sample-id ${prefix} \\
        ${args}

    # Compress and index the annotated VCF
    bgzip ${prefix}.vcf
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        htslib: \$(echo \$(htslib-config --version))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.vcf.gz
    touch ${prefix}_annotated.vcf.gz.tbi
    touch ${prefix}_annotated.tsv
    touch ${prefix}_annotation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        htslib: \$(echo \$(htslib-config --version))
    END_VERSIONS
    """
}
