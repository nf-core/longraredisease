// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test workflow
include { CALL_SNV } from '/home/nmahfel/nf-training/longraredisease/subworkflows/local/call_snv/tests/../main.nf'

workflow {

    // define custom rules for JSON that will be generated.
    def jsonOutput = createJsonOutput()
    def jsonWorkflowOutput = createJsonWorkflowOutput()

    def input = []

    // run dependencies
    

    // workflow mapping
    input = []
    
                input[0] = Channel.of([
                    [ id:'test_nanopore', platform:'ont' ],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/nanopore/bam/test.sorted.bam', checkIfExists: true),
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/nanopore/bam/test.sorted.bam.bai', checkIfExists: true)
                ])
                input[1] = Channel.of([
                    [ id:'reference' ],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta', checkIfExists: true)
                ])
                input[2] = Channel.of([
                    [ id:'reference' ],
                    file(params.modules_testdata_base_path + 'genomics/homo_sapiens/genome/genome.fasta.fai', checkIfExists: true)
                ])
                input[3] = false
                input[4] = Channel.empty()
                
    //----

    //run workflow
    CALL_SNV.run(input.toArray())
    
    if (CALL_SNV.output){

        // consumes all named output channels and stores items in a json file
        CALL_SNV.out.getNames().each { name ->
            serializeChannel(name, CALL_SNV.out.getProperty(name), jsonOutput, params.nf_test_output)
        }	  
    
        // consumes all unnamed output channels and stores items in a json file
        def array = CALL_SNV.out as List<Object>
        def i = 0
        array.each { output ->
            serializeChannel(i, output, jsonOutput, params.nf_test_output)
            i += 1
        }    	

    }

    workflow.onComplete = {

        def result = [
            success: workflow.success,
            exitStatus: workflow.exitStatus,
            errorMessage: workflow.errorMessage,
            errorReport: workflow.errorReport
        ]
        new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
        
    }
}


def serializeChannel(name, channel, jsonOutput, outputDir) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: { entry ->
            list.add(entry)
        },
        onComplete: {
            def map = new HashMap()
            map[_name] = list
            def filename = "${outputDir}/output_${_name}.json"
            new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}

def createJsonOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the process
    return [
        toJson: { obj ->
            def converted = convertPathsToStrings(obj)
            return groovy.json.JsonOutput.toJson(converted)
        }
    ]
}

def convertPathsToStrings(obj) {
    if (obj instanceof java.nio.file.Path) {
        return obj.toAbsolutePath().toString()
    } else if (obj instanceof Map) {
        return obj.collectEntries { k, v -> [k, convertPathsToStrings(v)] }
    } else if (obj instanceof Collection) {
        return obj.collect { it -> convertPathsToStrings(it) }
    } else {
        return obj
    }
}

def createJsonWorkflowOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the workflow
    return [
        toJson: { obj ->
            def filtered = removeNullValues(obj)
            return groovy.json.JsonOutput.toJson(filtered)
        }
    ]
}

def removeNullValues(obj) {
    if (obj instanceof Map) {
        return obj.findAll { _k, v -> v != null }.collectEntries { k, v -> [k, removeNullValues(v)] }
    } else if (obj instanceof Collection) {
        return obj.findAll { it -> it != null }.collect { it -> removeNullValues(it) }
    } else {
        return obj
    }
}