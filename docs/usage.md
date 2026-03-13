# nf-core/longraredisease: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/longraredisease/usage](https://nf-co.re/longraredisease/usage)

## Introduction

**nf-core/longraredisease** is a comprehensive Nextflow pipeline designed for rare disease diagnostics using Oxford Nanopore long-read sequencing data. The pipeline integrates multiple state-of-the-art tools for variant discovery, including:

- **Structural Variants (SVs):** Sniffles, CuteSV, SVIM with SURVIVOR merging
- **Single Nucleotide Variants (SNVs):** Clair3 and DeepVariant
- **Copy Number Variants (CNVs):** Spectre and HiFiCNV
- **Short Tandem Repeats (STRs):** Straglr
- **Methylation:** Modkit
- **Phasing:** LongPhase with haplotagging
- **Quality Control:** NanoPlot, mosdepth, and MultiQC

The pipeline supports singleton analysis as well as family-based (trio) analyses with variant annotation and phenotype-driven prioritization.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location:

```bash
--input '[path to samplesheet file]'
```

The samplesheet is a comma-separated file (CSV) with a header row and the following columns:

### Samplesheet format

| Column        | Required | Description                                                                          |
| ------------- | -------- | ------------------------------------------------------------------------------------ |
| `sample_id`   | Yes      | Unique sample identifier (no spaces)                                                 |
| `file_path`   | Yes      | Path to input file or directory (see details below)                                  |
| `hpo_terms`   | No       | HPO terms for phenotype annotation (format: `HP:0000001;HP:0000002`)                 |
| `sex`         | No       | Sex of the sample (`1`=male, `2`=female, `0`=unknown)                                |
| `phenotype`   | No       | Phenotype status (`1`=unaffected, `2`=affected, `0` or `-9`=missing)                 |
| `family_id`   | No       | Family identifier (required for trio analysis)                                       |
| `maternal_id` | No       | Sample ID of the maternal sample (must match another `sample_id` in the samplesheet) |
| `paternal_id` | No       | Sample ID of the paternal sample (must match another `sample_id` in the samplesheet) |

### Input file types

The `file_path` column can contain:

1. **Directory path**: A directory containing multiple files of the specified `input_type`
   - For `fastq` input: directory with `.fastq.gz` or `.fq.gz` files
   - For `ubam` input: directory with unaligned `.bam` files
2. **Single file path**: Path to a single file
   - For `fastq` input: single `.fastq.gz` or `.fq.gz` file
   - For `bam` or `ubam` input: single `.bam` file

The `input_type` is specified as a pipeline parameter (see [Workflow Options](#workflow-options)).

### Example samplesheets

#### Single sample (FASTQ input)

```csv title="samplesheet.csv"
sample_id,file_path,hpo_terms,sex,phenotype,family_id,maternal_id,paternal_id
sample1,/path/to/sample1.fastq.gz,HP:0002721;HP:0002110,1,2,,,
```

#### Trio analysis

```csv title="samplesheet_trio.csv"
sample_id,file_path,hpo_terms,sex,phenotype,family_id,maternal_id,paternal_id
proband,/path/to/proband.fastq.gz,HP:0002721;HP:0001263,1,2,family1,mother,father
mother,/path/to/mother.fastq.gz,,2,1,family1,,
father,/path/to/father.fastq.gz,,1,1,family1,,
```

#### Multiple families

```csv title="samplesheet_multi_family.csv"
sample_id,file_path,hpo_terms,sex,phenotype,family_id,maternal_id,paternal_id
fam1_child,/path/to/fam1_child/,HP:0001263,1,2,family1,fam1_mom,fam1_dad
fam1_mom,/path/to/fam1_mom/,,,family1,,
fam1_dad,/path/to/fam1_dad/,,,family1,,
fam2_child,/path/to/fam2_child.fastq.gz,HP:0002110,2,2,family2,fam2_mom,fam2_dad
fam2_mom,/path/to/fam2_mom.fastq.gz,,,family2,,
fam2_dad,/path/to/fam2_dad.fastq.gz,,,family2,,
```

Example samplesheets are provided in the `assets/` directory:

- [samplesheet_test_fastq.csv](../assets/samplesheet_test_fastq.csv)
- [samplesheet_test_ubam.csv](../assets/samplesheet_test_ubam.csv)

## Reference genome

The pipeline requires a reference genome in FASTA format. You can provide this using:

```bash
--fasta '/path/to/reference.fasta'
```

The pipeline supports common reference genomes (GRCh37/hg19, GRCh38/hg38). Ensure the reference genome is indexed (`.fai` file) or the pipeline will create the index automatically.

## Workflow options

### Input type

Specify the type of input files using the `--input_type` parameter:

```bash
--input_type 'fastq'   # FASTQ files (default: ubam)
--input_type 'ubam'    # Unaligned BAM files
--input_type 'bam'     # Aligned BAM files
```

### Sequencing platform

Specify the sequencing platform:

```bash
--sequencing_platform 'ont'     # Oxford Nanopore (default)
--sequencing_platform 'hifi'    # PacBio HiFi
--sequencing_platform 'pacbio'  # PacBio CLR
```

### Alignment options

By default, the pipeline uses Minimap2 for alignment. You can customize the alignment model or use Winnowmap instead:

```bash
--minimap2_model 'map-ont'     # Minimap2 preset (default: auto-detected)
--use_winnowmap true           # Use Winnowmap instead of Minimap2
--winnowmap_model 'map-ont'    # Winnowmap preset
--winnowmap_kmers '/path/to/repetitive_k15.txt'  # Repetitive k-mer file
```

### Analysis modules

Enable or disable specific analysis modules:

```bash
--skip_snv false             # Enable SNV calling (default: true)
--skip_sv false              # Enable SV calling (default: true)
--skip_cnv false             # Enable CNV calling (default: true)
--skip_str false             # Enable STR analysis (default: true)
--skip_methylation false     # Enable methylation calling (default: true)
--skip_phasing false         # Enable phasing (default: true)
```

### Trio analysis

Enable family-based analysis for samples with pedigree information:

```bash
--trio_analysis true         # Enable trio analysis (default: false)
--haplotag_bam true          # Haplotag BAM files with phase information (default: true)
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/longraredisease \
    --input ./samplesheet.csv \
    --outdir ./results \
    --fasta /path/to/reference.fasta \
    --input_type fastq \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

### Example commands

#### Basic singleton analysis

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results \
    --fasta GRCh38.fasta \
    --input_type fastq \
    -profile docker
```

#### Trio analysis with all modules enabled

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet_trio.csv \
    --outdir results_trio \
    --fasta GRCh38.fasta \
    --input_type ubam \
    --trio_analysis true \
    --haplotag_bam true \
    -profile docker
```

#### SV-focused analysis

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results_sv \
    --fasta GRCh38.fasta \
    --skip_snv true \
    --skip_cnv true \
    --skip_str true \
    --skip_methylation true \
    --skip_phasing true \
    -profile singularity
```

#### Analysis with target regions

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results_targeted \
    --fasta GRCh38.fasta \
    --filter_targets true \
    --targets_bed exome_targets.bed \
    -profile docker
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/nanoraredx -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/longraredisease
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/longraredisease releases page](https://github.com/nf-core/longraredisease/releases) and find the latest pipeline version - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
