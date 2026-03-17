process SPECTRE_CNVCALLER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/nourmahfel/spectre-ont:0.3.2"

    input:
    tuple val(meta), path(summary), path(regions_bed), path(regions_csi), path(vcf)
    tuple val(meta2), path(fasta)
    path(metadata_file)
    path(blacklist)
    val bin_size

    output:
    tuple val(meta), path("*.vcf")                       , emit: vcf
    tuple val(meta), path("*.bed")                       , emit: bed
    tuple val(meta), path("*.spc.gz")                    , emit: spc
    tuple val(meta), path("predicted_karyotype.txt")     , emit: txt
    tuple val(meta), path("windows_stats", type: 'dir')  , emit: winstats
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def coverage_dir = "mosdepth_${prefix}"

    """
    mkdir -p ${coverage_dir}

    # Stage the mosdepth files into the directory with proper naming
    cp -L ${summary} ${coverage_dir}/${meta.id}.mosdepth.summary.txt
    cp -L ${regions_bed} ${coverage_dir}/${meta.id}.regions.bed.gz
    cp -L ${regions_csi} ${coverage_dir}/${meta.id}.regions.bed.gz.csi

    # Now run spectre with the directory - use bin_size input parameter, not params
    spectre CNVCaller \\
        --coverage ${coverage_dir} \\
        --snv ${vcf} \\
        --reference ${fasta} \\
        --output-dir . \\
        --sample-id=${prefix} \\
        --metadata ${metadata_file} \\
        --blacklist ${blacklist} \\
        --bin-size ${bin_size} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spectre: \$(spectre --version 2>&1 | grep -oP 'version \\K[0-9.]+' || spectre --help 2>&1 | grep -oP 'v[0-9.]+' | sed 's/v//' || echo "1.0.0")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.bed
    touch ${prefix}.spc.gz
    touch predicted_karyotype.txt
    mkdir windows_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spectre: 0.3.2
    END_VERSIONS
    """
}
