process SNIFFLES_GENERATE_PLOTS {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/sniffles2-plot:0.2.1--54e822d08702535e"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}_plots", type: 'dir'), emit: plot_dir
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"

    """
    sniffles2_plot \\
        -i ${vcf} \\
        -o ${prefix}_plots \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles2_plot: \$(pip show sniffles2-plot 2>/dev/null | grep ^Version | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_plots
    touch ${prefix}_plots/plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles2_plot: \$(pip show sniffles2-plot 2>/dev/null | grep ^Version | sed 's/Version: //')
    END_VERSIONS
    """
}
