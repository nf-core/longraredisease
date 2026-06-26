process SNIFFLES_TRIO {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/bcftools_sniffles:1d48f8cb319ef79c'

    input:
        tuple val(meta), path(snf_files)  // Multiple SNF files as a list
        tuple val(meta2), path(reference)

    output:
        tuple val(meta), path("*.vcf"), emit: vcf
        path "versions.yml"           , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Join all SNF files with spaces for the command line
    def snf_input = snf_files.collect{ it.toString() }.join(' ')

    """
    sniffles \\
        --input ${snf_input} \\
        --vcf trio_tmp.vcf \\
        --reference $reference

    bcftools view -U -O v trio_tmp.vcf > "${prefix}.vcf"

    rm trio_tmp.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --version 2>&1 | head -n1 | sed 's/.*sniffles //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --version 2>&1 | head -n1 | sed 's/.*sniffles //')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
