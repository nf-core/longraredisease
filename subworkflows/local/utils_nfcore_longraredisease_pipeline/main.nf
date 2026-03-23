//
// Subworkflow with functionality specific to the nf-core/longraredisease pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )


    def command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

        def longraredisease_logo = """
    \033[0;34m _                                          _ _
    \033[0;34m| |                                        | (_)
    \033[0;34m| | ___  _ __   __ _ _ __ __ _ _ __ ___  __| |_ ___  ___  __ _ ___  ___
    \033[0;34m| |/ _ \\| '_ \\ / _` | '__/ _` | '__/ _ \\/ _` | / __|/ _ \\/ _` / __|/ _ \\
    \033[0;34m| | (_) | | | | (_| | | | (_| | | |  __/ (_| | \\__ \\  __/ (_| \\__ \\  __/
    \033[0;34m|_|\\___/|_| |_|\\__, |_|  \\__,_|_|  \\___|\\__,_|_|___/\\___|\\__,_|___/\\___|
    \033[0;34m                __/ |
    \033[0;34m               |___/
    \033[0;34m
    """

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        longraredisease_logo,
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    def workflowDependencies = [
        bam2fastq          : ["align"],
        qc                 : ["align"],
        bam_stats          : ["align"],
        mosdepth           : ["align"],
        multiqc_mosdepth   : ["mosdepth"],
        call_snv           : ["align"],
        annotate_snv       : ["call_snv"],
        haplotag_bam       : ["align", "call_snv"],
        call_sv            : ["align"],
        annotate_sv        : ["call_sv", "call_snv"],
        svanna_prioritise  : ["call_sv"],
        merge_sv           : ["call_sv"],
        trio_snv           : ["call_snv"],
        trio_sv            : ["call_sv"],
        call_str           : ["align"],
        call_cnv           : ["align", "mosdepth", "call_snv"],
        methyl             : ["align"],
        unify_vcf          : ["call_sv", "call_str", "call_cnv"],
        annotate_unified_vcf : ["unify_vcf"]
    ]

    def fileDependencies = [
        align            : ["fasta_file"],
        assembly         : ["fasta_file"],
        sambamba_depth   : ["sambamba_regions"],
        call_snv         : ["fasta_file"],
        snv_annotation   : ["snpeff_db"],
        call_sv          : ["fasta_file", "sniffles_tandem_file"],
        annotate_sv      : ["annotsv_annotations"],
        call_str         : ["straglr_bed"],
        str_annotation   : ["variant_catalogue"],
        call_cnv         : ["hificnv_exclude_bed", "hificnv_expected_cn_bed", "spectre_metadata", "spectre_blacklist"],
    ]

    def parameterStatus = [
        workflow: [
            // Map workflow flags to their negated parameter equivalents
            skip_call_snv        : !params.snv,
            skip_call_sv         : !params.sv,
            skip_methyl          : !params.methyl,
            skip_qc              : !params.qc,
            skip_call_str        : !params.str,
            skip_alignment       : (params.input_type == 'bam'),
        ],
        files: [
            fasta_file                  : params.fasta_file,
            sniffles_tandem_file        : params.sniffles_tandem_file,
            straglr_bed                 : params.straglr_bed,
            variant_catalogue           : params.variant_catalogue,
            hificnv_exclude_bed         : params.hificnv_exclude_bed,
            hificnv_expected_cn_bed     : params.hificnv_expected_cn_bed,
            spectre_metadata            : params.spectre_metadata,
            spectre_blacklist           : params.spectre_blacklist,
            annotsv_annotations         : params.annotsv_annotations,
            snpeff_db                   : params.snpeff_db,
            svanna_db                   : params.svanna_db,
            winnowmap_kmers             : params.winnowmap_kmers,
        ]
    ]

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(parameterStatus, workflowDependencies, fileDependencies)
    validateWorkflowCompatibility()

    def samplesheet_data = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")

    ch_samplesheet = Channel.fromList(samplesheet_data)
        .map { row ->
            // Handle the ArrayList structure from nf-schema
            if (row instanceof List) {
                def meta_map = row[0]
                def sample_id = meta_map.id ?: meta_map.toString()
                def meta = [
                    id: sample_id,
                    family_id: row[5] ?: 'unknown',
                    sex: row[3] ?: 0,
                    phenotype: row[4] ?: 0,
                    maternal_id: row[6] ?: "0",
                    paternal_id: row[7] ?: "0",
                    hpo_terms: row[2] ?: null
                ]
                def data = [
                    id: sample_id,
                    file_path: row[1],
                    hpo_terms: row[2] ?: null,
                    sex: row[3] ?: 0,
                    phenotype: row[4] ?: 0,
                    family_id: row[5] ?: 'unknown',
                    maternal_id: row[6] ?: "0",
                    paternal_id: row[7] ?: "0"
                ]
                return [meta, data]
            } else {
                error "Unexpected row type: ${row.getClass()}"
            }
        }

    validateWorkflowCompatibility()
    validateSVCallingParameters()

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  channel: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")


    def multiqc_reports = multiqc_report.ifEmpty([]).toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def validateInputParameters(statusMap, workflowDependencies, fileDependencies) {
    validateParameterCombinations(statusMap, workflowDependencies, fileDependencies)
}

//
// Validate channels from input samplesheet
//
def validateUniqueFilenamesPerSample(input) {
    // Filenames needs to be unique for each sample to avoid collisions when merging
    def fileNames = input[2].collect { input_path -> new File(input_path.toString()).name }
    if (fileNames.size() != fileNames.unique().size()) {
        error "Error: Input filenames needs to be unique for each sample."
    }
    return input
}

//
//
def validateUniqueSampleIDs(input) {
    def sample = input[0]
    def metas = input[1].collect()
    def families = metas.collect { meta -> meta.family_id }.unique()

    if (families.size() > 1) {
        error "Sample '${sample}' belongs to multiple families: ${families}. " +
                "Please make sure that there are no duplicate samples in the samplesheet."
    }
    return input
}

//
// Generate methods description for MultiQC
//
def methodsDescriptionText(mqc_methods_yaml) {
    // Convert to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def extractSoftwareFromVersions(module_yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def yamlData = yaml.load(module_yaml_file)
    // Extract all software (keys) from a module yaml

    def softwareInModule = yamlData.values().collect { software_and_version -> software_and_version.keySet() }.flatten()
    return softwareInModule
}

def extractSoftwareFromTopics(topics_channel) {
    topics_channel
        .map { toolBlockText ->
            toolBlockText
                .readLines()
                .drop(1) // Drop process name
                .collect { line -> line.trim().split(':')[0] }
        }
}

def generateReferenceHTML(tool_list, description) {
    def items = tool_list
        .collect { citation -> citation.trim() }
        .unique()                                // e.g. samtools and bcftools share citation
        .findAll { citation -> citation != "" }  // some tools does not have a citation, e.g. awk, gunzip

    if (description == 'citation') {
        return "  <p>Tools used in the workflow included: ${items.join(', ')}.</p>"
    } else if (description == 'bibliography') {
        return "  <h4>References</h4><ul><li>${items.join('</li><li>')}</li></ul>"
    }
}

def citationBibliographyText(ch_versions, ch_topic_versions_string, references_yaml, description) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def softwareReferences = yaml.load(references_yaml.text).tool

    def unwantedReferences = ['longraredisease', 'Nextflow']
    // These are not collected in ch_versions but should be referenced
    def baseTools = Channel.from(['nextflow', 'nf_core', 'bioconda', 'biocontainers', 'multiqc'])

    ch_versions
        .map { module_yaml -> extractSoftwareFromVersions(module_yaml) }
        .concat(extractSoftwareFromTopics(ch_topic_versions_string))
        .flatten() // split multi-tool modules
        .unique()
        .filter { tool -> !unwantedReferences.contains(tool) }
        .concat(baseTools)
        .collect { tool ->
            def toolDetails = softwareReferences[tool]
            if (toolDetails == null) {
                throw new IllegalStateException("Tool: '${tool}' not found in ${references_yaml}")
            }
            return toolDetails[description]
        }
        .sort()
        .map { tools -> generateReferenceHTML(tools, description) }
}


def validateParameterCombinations(statusMap, workflowDependencies, fileDependencies) {
    // Array to store errors
    def errors = []

    // For each of the "workflow", "files"
    statusMap.each { paramsType, paramsMap ->
        paramsMap.each { param, paramStatus ->
            if (paramsType == "files") {
                checkFileDependencies(param, fileDependencies, statusMap, errors)
            } else if (paramsType == "workflow") {
                checkWorkflowDependencies(param, workflowDependencies, statusMap, errors)
            }
        }
    }

    // Give error if there are any
    if (errors) {
        def error_string =
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  " + errors.join("\n  ") + "\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}


def checkWorkflowDependencies(String skip, Map workflowDependencies, Map statusMap, List errors) {
    // If the workflow is skipped, check if any dependent workflows are active
    def workflowIsSkipped = statusMap["workflow"][skip]
    if (!workflowIsSkipped) {
        return // Workflow is active, no error
    }

    // Find workflows that depend on this one
    def dependentWorkflows = workflowDependencies.findAll { workflow, dependencies ->
        dependencies.contains(skip.replace('skip_', '').replace('_calling', ''))
    }

    dependentWorkflows.each { workflow, dependencies ->
        def dependentSkip = "skip_${workflow}_calling"
        if (statusMap["workflow"][dependentSkip] == false) {
            errors << "Workflow '${workflow}' is active but depends on '${skip}' which is skipped"
        }
    }
}


def checkFileDependencies(String file, Map fileDependencies, Map statusMap, List errors) {
    def filePath = statusMap["files"][file]

    if (!filePath) {
        // Find workflows that require this file
        def requiringWorkflows = fileDependencies.findAll { workflow, files ->
            files.contains(file)
        }

        requiringWorkflows.each { workflow, files ->
            def workflowSkip = "skip_${workflow}_calling"
            def workflowIsActive = statusMap["workflow"][workflowSkip] == false

            if (workflowIsActive) {
                errors << "Workflow '${workflow}' is active but required file '${file}' is not provided"
            }
        }
    }
}

// Utility function to create channels from references
def createReferenceChannelFromPath(param, defaultValue = '', id = null) {
    return param ? Channel.fromPath(param, checkIfExists: true)
        .map { file_path -> [ [ id: id ?: file_path.simpleName ], file_path ] }
        .collect() : defaultValue
}

// Utility function to create channels from samplesheets
def createReferenceChannelFromSamplesheet(param, schema, defaultValue = '') {
    return param ? Channel.fromList(samplesheetToList(param, schema)) : defaultValue
}


def validateWorkflowCompatibility() {
    // Check CNV calling compatibility
    if (params.cnv && params.sequencing_platform in ['pacbio', 'hifi']) {
        if (!params.hificnv_exclude_bed || !params.hificnv_expected_cn_bed) {
            error "ERROR: HiFiCNV requires exclude bed and expected CN bed files. Please provide --hificnv_exclude_bed and --hificnv_expected_cn_bed parameters."
        }
    }

    // Check trio analysis requirements
    if (params.trio_analysis && (!params.snv && !params.sv)) {
        error "ERROR: Trio analysis requires either SNV or SV calling to be enabled."
    }
}

def validateSVCallingParameters() {
    // Only validate if SV calling is enabled
    if (!params.sv) {
        return
    }

    // Check if required SV parameters are provided
    if (params.sv && !params.sniffles_tandem_file) {
        error "ERROR: SV calling requires tandem repeat file. Please provide --sniffles_tandem_file parameter."
    }


}

def findKeysForValue(def valueToFind, Map map) {
    def keys = []

    map.each { entry ->
        def key = entry.key
        def value = entry.value

        if ((value instanceof List && value.contains(valueToFind)) || value == valueToFind) {
            keys << key
        }
    }
    return keys.isEmpty() ? null : keys
}

// Helper functions for relationship validation
def getParentalIds(samples, parental_id_type) {
    samples.collect { sample -> sample[parental_id_type] }.findAll { parental_id -> isNonZeroNonEmpty(parental_id) }
}

def addRelationshipsToMeta(samples) {
    // This function adds relationships to the samples based on their parental IDs.
    def maternal_ids = getParentalIds(samples, 'maternal_id')
    def paternal_ids = getParentalIds(samples, 'paternal_id')
    def parents_ids = maternal_ids + paternal_ids
    def grandparents_ids = samples.findAll { sample -> sample.id in parents_ids }.collect { sample -> sample.maternal_id } +
                            samples.findAll { sample -> sample.id in parents_ids }.collect { sample -> sample.paternal_id }

    samples.each { sample ->
        sample.relationship = sample.id in grandparents_ids ? 'unknown' :
                                sample.id in maternal_ids ? 'mother' :
                                sample.id in paternal_ids ? 'father' :
                                isChild(sample, maternal_ids, paternal_ids) ? 'child' : 'unknown'

        sample.two_parents = isChildWithTwoParents(sample, maternal_ids, paternal_ids)

        // Find children of this specific parent
        sample.children = []
        sample.has_other_parent = false

        if (isParent(sample)) {
            def children = getChildrenForParent(samples, sample.id)
            sample.children = children.collect { meta -> meta.id }

            if (isMother(sample)) {
                sample.has_other_parent = children.any { child -> hasFather(child, paternal_ids) }
            } else if (isFather(sample)) {
                sample.has_other_parent = children.any { child -> hasMother(child, maternal_ids) }
            }
        }
    }
    return samples
}

def getChildrenForParent(samples, parent_id) {
    samples.findAll { sample -> sample.maternal_id == parent_id || sample.paternal_id == parent_id }
}

def isChild(sample, maternal_ids, paternal_ids) {
    hasMother(sample, maternal_ids) || hasFather(sample, paternal_ids)
}

def isChildWithTwoParents(sample, maternal_ids, paternal_ids) {
    hasMother(sample, maternal_ids) && hasFather(sample, paternal_ids)
}

def hasMother(sample, maternal_ids) {
    sample.maternal_id in maternal_ids
}

def hasFather(sample, paternal_ids) {
    sample.paternal_id in paternal_ids
}

def isFemale(sample) {
    sample.sex == 2
}

def isMale(sample) {
    sample.sex == 1
}

def isMother(sample) {
    sample.relationship == 'mother'
}

def isFather(sample) {
    sample.relationship == 'father'
}

def isParent(sample) {
    isMother(sample) || isFather(sample)
}

def boolean isNonZeroNonEmpty(value) {
    (value instanceof String && value != "" && value != "0") ||
    (value instanceof Number && value != 0)
}
