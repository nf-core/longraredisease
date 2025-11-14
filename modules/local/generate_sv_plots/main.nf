process SV_PLOT {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/nourmahfel1/sniffles2_plot"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}"), emit: plot_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_dir = "${prefix}"

    """
    # Set matplotlib config directory to a writable location
    export MPLCONFIGDIR=\${PWD}/matplotlib_config
    mkdir -p \${MPLCONFIGDIR}

    python3 -m sniffles2_plot \\
        -i ${vcf} \\
        -o ${output_dir} \\
        ${args}

    # Ensure output directory has content
    touch ${output_dir}/.keep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles2_plot: \$(python3 -m sniffles2_plot --version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "unknown")
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_dir = "${prefix}_sv_plot_output"

    """
    mkdir -p ${output_dir}
    touch ${output_dir}/plot1.png
    touch ${output_dir}/plot2.png
    touch ${output_dir}/summary.txt
    touch ${output_dir}/.keep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles2_plot: unknown
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
