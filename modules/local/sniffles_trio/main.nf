// process SNIFFLES {
//     tag "$meta.id"
//     label 'process_high'

//     conda "${moduleDir}/environment.yml"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/sniffles:2.4--pyhdfd78af_0' :
//         'biocontainers/sniffles:2.4--pyhdfd78af_0' }"

//     input:
//     tuple val(meta), path(input), path(index)
//     tuple val(meta2), path(fasta)
//     tuple val(meta3), path(tandem_file)
//     val(vcf_output)
//     val(snf_output)


//     output:
//     tuple val(meta), path("*.vcf.gz")    , emit: vcf, optional: true
//     tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
//     tuple val(meta), path("*.snf")       , emit: snf, optional: true
//     path "versions.yml"                  , emit: versions

//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     def args = task.ext.args ?: ''
//     def prefix = task.ext.prefix ?: "${meta.id}"
//     def reference = fasta ? "--reference ${fasta}" : ""
//     def tandem_repeats = tandem_file ? "--tandem-repeats ${tandem_file}" : ''
//     def vcf = vcf_output ? "--vcf ${prefix}.vcf.gz": ''
//     def snf = snf_output ? "--snf ${prefix}.snf": ''

//     """
//     sniffles \\
//         --input $input \\
//         $reference \\
//         -t $task.cpus \\
//         $tandem_repeats \\
//         $vcf \\
//         $snf \\
//         $args

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         sniffles: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version //')
//     END_VERSIONS
//     """

//     stub:
//     def prefix = task.ext.prefix ?: "${meta.id}"
//     def vcf = vcf_output ? "echo \"\" | gzip > ${prefix}.vcf.gz; touch ${prefix}.vcf.gz.tbi": ''
//     def snf = snf_output ? "touch ${prefix}.snf": ''

//     """
//     ${vcf}
//     ${snf}

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         sniffles: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version //')
//     END_VERSIONS
//     """
// }
// // Trio SV

// // Multi-sample SV calling
// process sniffles2_joint {
//     cpus 4
//     memory 6.GB
//     label "wf_human_sv"
//     input:
//         tuple val(joint_prefix), path("snfs/*")
//         tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
//         val(suffix) // Extra label to add to the output files
//     output:
//         tuple val(joint_meta), path("${joint_prefix}.${suffix}.vcf"), emit: vcf
//     script:
//         // Perform internal phasing only if snp not requested; otherwise, use joint phasing.
//         def phase = params.phased ? "--phase" : ""
//         def min_sv_len = params.min_sv_length ? "--minsvlen ${params.min_sv_length}" : ""
//         def sniffles_args = params.sniffles_args ?: ''
//         joint_meta = ['family_id': joint_prefix]
//     """
//     sniffles \
//         --threads $task.cpus \
//           ${min_sv_len} \
//         --input snfs/* \
//         --vcf "${joint_prefix}.trio.sv.tmp.vcf" \
//         --reference $ref \
//         $phase \
//         $sniffles_args
//     bcftools view -U -O v "${joint_prefix}.trio.sv.tmp.vcf" > "${joint_prefix}.${suffix}.vcf"
//     """
// }

// // TO DO: Update wf-hum-var sniffles2 to accept tuple for multisample
// // NOTE VCF entries for alleles with no support are removed to prevent them from
// //      breaking downstream parsers that do not expect them
// // --input-exclude-flags 2308: Remove unmapped (4), non-primary (256) and supplemental (2048) alignments
// process sniffles2 {
//     label "wf_human_sv"
//     cpus 4
//     memory 24.GB
//     input:
//         tuple path(xam), path(xam_idx), val(xam_meta),
//             path(tr_bed),
//             path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
//         val suffix // Extra label to add to the output files
//     output:
//         tuple val(xam_meta), path("${xam_meta.alias}.${suffix}.tmp.vcf"), emit: vcf
//         tuple val(xam_meta), path("${xam_meta.alias}.${suffix}.vcf.gz"), path("${xam_meta.alias}.${suffix}.vcf.gz.tbi"), emit: compressed
//         tuple val(xam_meta), path("${xam_meta.alias}.${suffix}.snf"), emit: snf
//     script:
//         def tr_arg = ""
//         if (tr_bed.name != 'OPTIONAL_FILE'){
//             tr_arg = "--tandem-repeats ${tr_bed}"
//         } else {
//             log.info "Automatically selecting TR BED: hg38.trf.bed"
//             tr_arg = "--tandem-repeats \${WFSV_TRBED_PATH}/hg38.trf.bed"
//         }
//         def sniffles_args = params.sniffles_args ?: ''
//         def min_sv_len = params.min_sv_length ? "--minsvlen ${params.min_sv_length}" : ""
//         // Perform internal phasing only if snp not requested; otherwise, use joint phasing.
//         def phase = params.phased ? "--phase" : ""
//     """
//     sniffles \
//         --threads $task.cpus \
//         --sample-id ${xam_meta.alias} \
//         --output-rnames \
//         ${min_sv_len} \
//         --cluster-merge-pos $params.cluster_merge_pos \
//         --input $xam \
//         --reference $ref \
//         --input-exclude-flags 2308 \
//         --snf "${xam_meta.alias}.${suffix}.snf" \
//         $tr_arg \
//         $sniffles_args \
//         $phase \
//         --vcf "${xam_meta.alias}.sniffles.vcf"
//     # After running sniffles filter out uncalled and bgzip and index for individual vcfs
//     # unfiltered version required for refine snp with sv process
//     bcftools view --exclude-uncalled --output-type v "${xam_meta.alias}.sniffles.vcf" > "${xam_meta.alias}.${suffix}.tmp.vcf"
//     cp "${xam_meta.alias}.${suffix}.tmp.vcf" "${xam_meta.alias}.${suffix}.vcf"
//     bgzip "${xam_meta.alias}.${suffix}.vcf"
//     tabix --force --preset vcf "${xam_meta.alias}.${suffix}.vcf.gz"
//     """
// }


// // TO DO: Update this in humvar to accept suffix
// // NOTE This is the last touch the VCF has as part of the workflow,
// //  we'll rename it with its desired output name here
// process sortVCF {
//     label "wf_human_sv"
//     cpus 2
//     memory 4.GB
//     input:
//         tuple val(xam_meta), path(vcf)
//         val suffix
//     output:
//         tuple val(xam_meta), path("${xam_meta.family_id}.${suffix}.vcf.gz"), emit: vcf_gz
//         tuple val(xam_meta), path("${xam_meta.family_id}.${suffix}.vcf.gz.tbi"), emit: vcf_tbi
//     script:
//     """
//     bcftools sort -m 2G -T ./ -O z $vcf > "${xam_meta.family_id}.${suffix}.vcf.gz"
//     tabix -p vcf "${xam_meta.family_id}.${suffix}.vcf.gz"
//     """
// }

// // Run on individual VCF and joint VCF as that is made from SNFS which have not been filtered
// process filterCalls {
//     label "wf_human_sv"
//     cpus 4
//     memory 4.GB
//     input:
//         tuple val(xam_meta), path("filter_calls.vcf.gz"), path("filter_calls.vcf.gz.tbi"), path("input.bed")
//         val(chromosome_codes)
//         val suffix
//     output:
//         tuple val(xam_meta), path("${xam_meta.alias}.${suffix}.vcf.gz"), path("${xam_meta.alias}.${suffix}.vcf.gz.tbi"), emit: vcf
//     script:
//     String ctgs = chromosome_codes.join(',')
//     def ctgs_filter = params.include_all_ctgs ? "" : "-r ${ctgs}"
//     String bed = params.bed ? "-T input.bed  --targets-overlap 1" : ""
//     """
//     bcftools view --threads ${task.cpus} ${ctgs_filter} ${bed} "filter_calls.vcf.gz" > "${xam_meta.alias}.${suffix}.vcf"
//     bgzip "${xam_meta.alias}.${suffix}.vcf"
//     tabix -p vcf "${xam_meta.alias}.${suffix}.vcf.gz"
//     """
// }

// //TO DO: Make wf-hum-var report process and python script
// //reusable currently report name is hard coded
// process report {
//     label "wf_common"
//     publishDir "${params.out_dir}/${xam_meta.alias}", mode: 'copy', pattern: "*wf-trio-sv-report.html"
//     cpus 1
//     memory 4.GB
//     input:
//         tuple val(xam_meta), path(vcf), path(tbi)
//         file versions
//         path "params.json"
//     output:
//         path "${xam_meta.alias}.wf-trio-sv-report.html", emit: html
//     script:
//         String workflow_name = workflow.manifest.name.replace("epi2me-labs/", "")
//         def report_name = "${xam_meta.alias}.wf-trio-sv-report.html"
//     """

//     workflow-glue report_sv \
//         $report_name \
//         --workflow_name ${workflow_name} \
//         --vcf $vcf \
//         --params params.json \
//         --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
//         --versions $versions \
//         --revision ${workflow.revision} \
//         --commit ${workflow.commitId} \
//         --output_json "${xam_meta.alias}.svs.json" \
//         --workflow_version ${workflow.manifest.version}
//     """
// }


// process getVersions {
//     label "wf_human_sv"
//     cpus 1
//     memory 2.GB
//     output:
//         path "versions.txt"
//     script:
//     """
//     trap '' PIPE # suppress SIGPIPE without interfering with pipefail
//     sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt
//     bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
//     samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
//     """
// }



// process makeJointReport {
//     label "wf_common"
//     publishDir "${params.out_dir}", mode: 'copy', pattern: "*.wf-trio-sv-report.html"
//     cpus 1
//     memory 4.GB
//     input:
//         tuple val(family_id), path(rtg_mendelian)
//         path versions
//         path "params.json"
//         path "ped_file.ped"
//     output:
//         path "${family_id}.wf-trio-sv-report.html", emit: 'report'
//     script:
//         def report_name = "${family_id}.wf-trio-sv-report.html"
//         def wfversion = workflow.manifest.version
//         if( workflow.commitId ){
//             wfversion = workflow.commitId
//         }
//         """
//         workflow-glue report_joint_sv \
//         $report_name \
//         --versions $versions \
//         --params params.json \
//         --sample_name $family_id \
//         --wf_version ${workflow.manifest.version} \
//         --rtg_mendelian ${rtg_mendelian} \
//         --ped_file "ped_file.ped"
//         """
// }
