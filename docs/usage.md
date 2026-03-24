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

nextflow run nf-core/longraredisease \
 --input ./samplesheet.csv \
 --outdir ./results \
 --fasta /path/to/reference.fasta \
 --input_type fastq \
 -profile docker

````

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
````

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

nextflow run nf-core/nanoraredx -profile docker -params-file params.yaml

````

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
````

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/longraredisease
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/longraredisease releases page](https://github.com/nf-core/longraredisease/releases) and find the latest pipeline version - numeric only (eg. `1.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0.0`. Of course, you can switch to another version by changing the number after the `-r` flag.
