process DOWNSAMPLE_SV {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/bcftools_pip_confargparse:4f3c18aa8341a070"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(mosdepth_summary)
    tuple val(meta3), path(target_bed)
    val chromosome_codes
    val min_read_support
    val min_read_support_limit

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: filterbycov_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_name = "${prefix}"
    def ctgs = chromosome_codes.join(',')
    def ctgs_filter = "--contigs ${ctgs}"

    """
    # Filter bed file for callable regions - handle gzipped files directly
    # Check if target_bed is provided and is a real file (not empty/null)
    if [[ -f "${target_bed}" && "${target_bed}" != "OPTIONAL_FILE" ]]; then
        echo "Processing BED file: ${target_bed}" >&2

        if [[ "${target_bed}" == *.gz ]]; then
            zcat ${target_bed} | awk '\$4 == "CALLABLE" || \$4 == "HIGH_COVERAGE"' > callable_regions.bed
        else
            awk '\$4 == "CALLABLE" || \$4 == "HIGH_COVERAGE"' ${target_bed} > callable_regions.bed
        fi

        # Check if the filtered BED file has any content
        if [[ -s callable_regions.bed ]]; then
            echo "Found \$(wc -l < callable_regions.bed) callable regions" >&2
            target_bed_arg="--target_bedfile callable_regions.bed"
        else
            echo "Warning: No CALLABLE or HIGH_COVERAGE regions found in BED file. Skipping region filtering." >&2
            target_bed_arg=""
        fi
    else
        echo "No target BED file provided. Skipping region filtering." >&2
        target_bed_arg=""
    fi

    # Use input VCF directly (already compressed with index)
    input_vcf="${vcf}"

    # Verify index exists and is accessible
    if [[ ! -f "${tbi}" ]]; then
        echo "Error: Index file ${tbi} not found"
        exit 1
    fi

    # Extract average depth from mosdepth summary
    AVG_DEPTH=\$(awk '\$1 == "total" {print \$4}' ${mosdepth_summary})
    echo "Average depth: \$AVG_DEPTH" >&2

    # Generate filtering command with PASS filter option
    epi2me_downsample_sv.py \\
        --bcftools_threads ${task.cpus} \\
        \$target_bed_arg \\
        --vcf \$input_vcf \\
        --depth_summary ${mosdepth_summary} \\
        --min_read_support ${min_read_support} \\
        --min_read_support_limit ${min_read_support_limit} \\
        ${ctgs_filter} \\
        ${args} > filter_command.sh

    # Show the generated command for debugging
    echo "Generated filter command:" >&2
    cat filter_command.sh >&2

    # Execute filtering and compress output with custom naming
    bash filter_command.sh | bcftools view -O z -o ${output_name}.vcf.gz

    # Index the output VCF
    tabix -p vcf ${output_name}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //g')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "covFiltered"
    def output_name = "${prefix}_${suffix}"
    """
    touch ${output_name}.vcf.gz
    touch ${output_name}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //g')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
