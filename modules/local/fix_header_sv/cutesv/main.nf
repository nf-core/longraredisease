process RE2SUPPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.8"
    container "community.wave.seqera.io/library/sed:4.9--b22139a895c82f4b"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}_normalized.vcf"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Rename INFO header ID=RE -> ID=SUPPORT
    # Rename INFO key occurrences RE= -> SUPPORT= for all records (non-header lines)
    sed -e 's/^##INFO=<ID=RE,/##INFO=<ID=SUPPORT,/' \\
        -e '/^#/!s/\\(^\\|;\\)RE=/\\1SUPPORT=/g' \\
        "${vcf}" > "${prefix}_normalized.vcf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | head -n1 | sed 's/^.*sed //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_normalized.vcf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | head -n1 | sed 's/^.*sed //; s/ .*\$//')
    END_VERSIONS
    """
}
