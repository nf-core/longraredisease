process NORMALIZE_JASMINE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gawk=5.1.0"
    container "community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.vcf"), path("${prefix}_samples.txt"), emit: vcf_with_samples
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Input validation
    [ -f "${vcf}" ] || { echo "Input not found: ${vcf}" >&2; exit 1; }

    # Step 1: Replace any REF '.' with 'N'
    if [[ "${vcf}" == *.gz ]]; then
        zcat "${vcf}" | awk 'BEGIN{FS=OFS="\\t"}
            /^#/ {print; next}
            \$4=="." { \$4="N" }
            {print}
        ' > __tmp.refN.vcf
    else
        awk 'BEGIN{FS=OFS="\\t"}
            /^#/ {print; next}
            \$4=="." { \$4="N" }
            {print}
        ' "${vcf}" > __tmp.refN.vcf
    fi

    # Step 2: Create header definitions file
    cat > __add_headers.txt <<'HDR'
##ALT=<ID=TRA,Description="Translocation event (symbolic allele)">
##ALT=<ID=BND,Description="Breakend (symbolic allele)">
##ALT=<ID=INV,Description="Inversion (symbolic allele)">
##INFO=<ID=STRANDS,Number=1,Type=String,Description="Strand orientation for adjacencies (e.g., ++,+-,-+,--)">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Read support">
##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Long-read support">
##INFO=<ID=SUPPORT_SA,Number=1,Type=Integer,Description="Number of supplementary/secondary alignments supporting the variant">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">
##INFO=<ID=COVERAGE,Number=1,Type=Integer,Description="Coverage at site">
##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Std dev of variant length across callers">
##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Std dev of breakpoint position across callers">
##INFO=<ID=AVG_VARCALLS,Number=1,Type=Float,Description="Average number of variant calls">
##INFO=<ID=AVG_LEN,Number=1,Type=Float,Description="Average length">
##INFO=<ID=AVG_START,Number=1,Type=Float,Description="Average start position">
##INFO=<ID=AVG_END,Number=1,Type=Float,Description="Average end position">
##INFO=<ID=STARTVARIANCE,Number=1,Type=Float,Description="Start position variance">
##INFO=<ID=ENDVARIANCE,Number=1,Type=Float,Description="End position variance">
##INFO=<ID=VARCALLS,Number=1,Type=Integer,Description="Number of variant calls">
##INFO=<ID=ALLVARS_EXT,Number=.,Type=String,Description="All variants external">
##INFO=<ID=SUPP_VEC_EXT,Number=1,Type=String,Description="Support vector external">
##INFO=<ID=IDLIST_EXT,Number=.,Type=String,Description="ID list external">
##INFO=<ID=SUPP_EXT,Number=1,Type=Integer,Description="Support external">
##INFO=<ID=INTRASAMPLE_IDLIST,Number=.,Type=String,Description="Intrasample ID list">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for second breakpoint">
##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description="Genomic origin of interspersed duplication seems to be deleted">
##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description="Standard deviation in span of merged SV signatures">
##INFO=<ID=STD_POS,Number=1,Type=Float,Description="Standard deviation in position of merged SV signatures">
##INFO=<ID=STD_POS1,Number=1,Type=Float,Description="Standard deviation of breakend 1 position">
##INFO=<ID=STD_POS2,Number=1,Type=Float,Description="Standard deviation of breakend 2 position">
##INFO=<ID=SEQS,Number=.,Type=String,Description="Insertion sequences from all supporting reads">
HDR

    # Step 3: Add headers before #CHROM line
    awk -v addfile="__add_headers.txt" '
        BEGIN{
        while ((getline l < addfile) > 0) need[l]=1
        }
        /^##/ { hdr[\$0]=1; print; next }
        /^#CHROM/ {
        for (l in need) if (!(l in hdr)) print l
        print; next
        }
        { print }
    ' __tmp.refN.vcf > __tmp.hdr.vcf

    # Step 4: Remove END field from BND variants
    awk 'BEGIN{FS=OFS="\\t"}
        /^#/ {print; next}
        {
        info=\$8
        if (info ~ /(^|;)SVTYPE=BND(\$|;)/) {
            gsub(/(^|;)END=[^;]*/, "", info)
            gsub(/;;+/, ";", info)
            sub(/^;/, "", info); sub(/;\$/, "", info)
            if (info == "") info="."
            \$8=info
        }
        print
        }' __tmp.hdr.vcf > "${prefix}.vcf"

    # Step 5: Create sample file for reheader
    echo "${meta.id}" > "${prefix}_samples.txt"

    # Clean up
    rm -f __tmp.refN.vcf __tmp.hdr.vcf __add_headers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version 2>&1 | head -n1 | sed 's/^.*gawk //; s/,.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.vcf"
    echo "${meta.id}" > "${prefix}_samples.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version 2>&1 | head -n1 | sed 's/^.*gawk //; s/,.*\$//')
    END_VERSIONS
    """
}
