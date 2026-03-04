process SNIFFLES_TRIO {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/bcftools_sniffles:1d48f8cb319ef79c'

    input:
        tuple val(meta), path(snf_files)  // Multiple SNF files as a list
        tuple val(meta2), path(reference)

    output:
        tuple val(meta), path("*.vcf"), emit: vcf

    script:
    def phase = params.phased ? "--phase" : ""
    def sniffles_args = params.sniffles_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Join all SNF files with spaces for the command line
    def snf_input = snf_files.collect{ it.toString() }.join(' ')

    """
    sniffles \\
        --input ${snf_input} \\
        --vcf trio_tmp.vcf \\
        --reference $reference \\
        $phase \\
        $sniffles_args

    bcftools view -U -O v trio_tmp.vcf > "${prefix}.vcf"

    rm trio_tmp.vcf
    """
}
