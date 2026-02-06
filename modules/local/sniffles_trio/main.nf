process sniffles2_joint {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/bcftools_sniffles:1d48f8cb319ef79c'

    input:
        tuple val(meta), path("snfs/*")
        tuple val(meta2), path(reference)

    output:
        tuple val(meta), path("*.vcf"), emit: vcf

    script:
    def phase = params.phased ? "--phase" : ""
    def sniffles_args = params.sniffles_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sniffles \
        --input snfs/* \
        --vcf "${prefix}.trio.sv.tmp.vcf" \
        --reference $reference \
        $phase \
        $sniffles_args
    bcftools view -U -O v "${prefix}.trio.sv.tmp.vcf" > "${prefix}.vcf"
    """
}
