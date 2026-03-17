<h1>
    <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-longraredisease_logo_dark.png">
    <img alt="nf-core/longraredisease" src="docs/images/nf-core-longraredisease_logo_light.png">
    </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/longraredisease/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/longraredisease/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/longraredisease/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/longraredisease/actions/workflows/linting.yml)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/longraredisease/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/longraredisease)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23longraredisease-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/longraredisease)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

---

## Introduction

**nf-core/longraredisease** is a specialized bioinformatics pipeline for **structural variant (SV) detection and clinical interpretation** from long-read sequencing data (Oxford Nanopore and PacBio). Designed for rare disease diagnostics, it delivers high-confidence variant discovery through multi-caller consensus, family-based analysis, and phenotype-driven prioritization.

### 🎯 **Primary Focus: Structural Variant Detection**

The pipeline excels at identifying and interpreting structural variants through:

- **Multi-caller SV consensus** - Sniffles, CuteSV, SVIM with JASMINE merging
- **Phase-aware calling** - Haplotype-resolved SV detection using LongPhase
- **Family analysis** - Trio-based joint calling and de novo variant detection
- **Clinical annotation** - AnnotSV with disease database integration
- **Phenotype prioritization** - SVANNA-based ranking using HPO terms

### 📊 **Analysis Capabilities**

**Core SV Analysis (Always Enabled):**

- ✅ **Structural Variants** - Multi-caller detection (DEL, INS, DUP, INV, BND)
- ✅ **Phasing** - Long-range haplotyping with LongPhase
- ✅ **Quality Control** - Comprehensive QC with NanoPlot, mosdepth, MultiQC

**Optional Analyses:**

- 🧬 **Single Nucleotide Variants** - Clair3 or DeepVariant (enable with `--snv true`)
- 📈 **Copy Number Variants** - Spectre or HiFiCNV (enable with `--cnv true`)
- 🔁 **Short Tandem Repeats** - Straglr genotyping (enable with `--str true`)
- 🧪 **DNA Methylation** - Modkit extraction for ONT (enable with `--methyl true`)

---

## Requirements

### Software

- **Nextflow** ≥25.04.0 (DSL2)
- **Container engine:** Docker, Singularity/Apptainer, or Podman
- **Java** ≥17 (required by Nextflow)

### Recommended Hardware

| Analysis Type         | CPU Cores | Memory   | Storage |
| --------------------- | --------- | -------- | ------- |
| **Single WGS sample** | 8-16      | 32-64 GB | 100 GB  |

**Notes:**

- Coverage recommendations: ≥10x for accurate SV calling, ≥30x for high-confidence trio analysis
- Storage includes space for input data, intermediate files, and results
- Adjust `--max_cpus` and `--max_memory` parameters based on available resources

---

## Quick Start

### 1. Install Nextflow

```bash
# Install Nextflow (≥25.04.0)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

### 2. Test the Pipeline

```bash
# Run with test data
nextflow run nf-core/longraredisease \
    -profile test,docker \
    --outdir test_results
```

### 3. Run with the Longraredisease Test Data

**Minimal SV-focused run:**

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results \
    --fasta reference.fasta \
    --sequencing_platform ont \
    -profile docker
```

**With family analysis and phenotype prioritization:**

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results \
    --fasta reference.fasta \
    --sequencing_platform ont \
    --trio_analysis true \
    --run_svanna true \
    --svanna_db /path/to/svanna_db \
    -profile docker
```

See [docs/usage.md](docs/usage.md) for complete examples and parameter details.

---

## Input Requirements

### Required Inputs

| Parameter               | Description                      | Format                      | Example           |
| ----------------------- | -------------------------------- | --------------------------- | ----------------- |
| `--input`               | Samplesheet with sample metadata | CSV                         | `samplesheet.csv` |
| `--outdir`              | Output directory                 | Path                        | `./results`       |
| `--fasta`               | Reference genome FASTA           | `.fasta`/`.fa`              | `GRCh38.fasta`    |
| `--sequencing_platform` | Platform type                    | `ont` or `pacbio` or `hifi` | `ont`             |

### Samplesheet Format

The input samplesheet is a CSV file with the following columns:

**Minimal format (single samples):**

```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,/path/to/sample2.bam.bai
```

**Family analysis format (trios):**

```csv
sample,bam,bai,family,paternal_id,maternal_id,sex,phenotype,hpo_terms
proband,proband.bam,proband.bam.bai,family1,father,mother,1,affected,"HP:0001250,HP:0002066"
father,father.bam,father.bam.bai,family1,0,0,1,unaffected,
mother,mother.bam,mother.bam.bai,family1,0,0,2,unaffected,
```

**Column descriptions:**

- `sample` - Unique sample identifier
- `bam` - Path to aligned BAM file
- `bai` - Path to BAM index file
- `family` - Family identifier (for trio analysis)
- `paternal_id` - Father's sample ID (or `0` if not in study)
- `maternal_id` - Mother's sample ID (or `0` if not in study)
- `sex` - `1` = male, `2` = female, `0` = unknown
- `phenotype` - `affected` or `unaffected`
- `hpo_terms` - Comma-separated HPO terms (e.g., `HP:0001250,HP:0002066`)

### Optional Inputs

| Parameter      | Description             | Required For             |
| -------------- | ----------------------- | ------------------------ |
| `--bed`        | Target regions BED file | Targeted sequencing      |
| `--annotsv_db` | AnnotSV database path   | SV annotation            |
| `--svanna_db`  | SVANNA database path    | Phenotype prioritization |
| `--str_bed`    | STR loci BED file       | STR analysis             |

---

## Key Parameters

### Core Analysis Toggles

**Structural variant analysis is always enabled.** Optional analyses:

| Parameter  | Description                                   | Default |
| ---------- | --------------------------------------------- | ------- |
| `--snv`    | Enable SNV calling (Clair3/DeepVariant)       | `false` |
| `--cnv`    | Enable CNV detection (Spectre)                | `false` |
| `--str`    | Enable STR genotyping (Straglr)               | `false` |
| `--methyl` | Enable methylation calling (Modkit, ONT only) | `false` |

### SV Detection Parameters

| Parameter            | Description                               | Default |
| -------------------- | ----------------------------------------- | ------- |
| `--run_cutesv`       | Enable CuteSV caller                      | `true`  |
| `--run_svim`         | Enable SVIM caller (recommended for BNDs) | `false` |
| `--haplotag_bam`     | Haplotag BAM for phase-aware SV calling   | `true`  |
| `--min_sv_size`      | Minimum SV size to report (bp)            | `30`    |
| `--min_read_support` | Minimum supporting reads                  | `auto`  |

### Family Analysis Parameters

| Parameter         | Description                            | Default |
| ----------------- | -------------------------------------- | ------- |
| `--trio_analysis` | Enable trio/family-based calling       | `false` |
| `--run_svanna`    | Enable phenotype-driven prioritization | `false` |
| `--svanna_db`     | Path to SVANNA database                | -       |

### Multi-caller Consensus Parameters

| Parameter               | Description                               | Default |
| ----------------------- | ----------------------------------------- | ------- |
| `--jasmine_max_dist`    | Max distance for merging breakpoints (bp) | `1000`  |
| `--jasmine_min_support` | Min callers supporting merged variant     | `2`     |
| `--jasmine_spec_reads`  | Min supporting reads for consensus        | `3`     |

### Platform-specific Settings

| Parameter               | Description               | Options                         |
| ----------------------- | ------------------------- | ------------------------------- |
| `--sequencing_platform` | Sequencing platform       | `ont`, `pacbio`                 |
| `--preset`              | Minimap2 alignment preset | `map-ont`, `map-hifi`, `map-pb` |
| `--snv_caller`          | SNV caller choice         | `clair3`, `deepvariant`         |

---

## Usage Examples

### 1. Standard SV Analysis (Single Sample)

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results \
    --fasta GRCh38.fasta \
    --sequencing_platform ont \
    -profile docker
```

### 2. Comprehensive Analysis (SVs + SNVs + CNVs)

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir results \
    --fasta GRCh38.fasta \
    --sequencing_platform pacbio \
    --snv true \
    --cnv true \
    --str true \
    -profile singularity
```

### 3. Family Trio Analysis with Phenotype Prioritization

```bash
nextflow run nf-core/longraredisease \
    --input trio_samplesheet.csv \
    --outdir family_results \
    --fasta GRCh38.fasta \
    --sequencing_platform ont \
    --trio_analysis true \
    --run_svanna true \
    --svanna_db /databases/svanna_data \
    --annotsv_db /databases/AnnotSV \
    -profile docker
```

### 4. High-Sensitivity SV Detection

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir sensitive_results \
    --fasta GRCh38.fasta \
    --sequencing_platform ont \
    --run_svim true \
    --min_sv_size 20 \
    --min_read_support 2 \
    --jasmine_min_support 1 \
    -profile docker
```

### 5. Targeted Sequencing with BED File

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --outdir targeted_results \
    --fasta GRCh38.fasta \
    --bed targets.bed \
    --sequencing_platform ont \
    -profile docker
```

---

## Output Structure

```
results/
├── pipeline_info/              # Pipeline execution reports
│   ├── execution_report.html   # Resource usage timeline
│   ├── execution_timeline.html # Process execution graph
│   └── multiqc_report.html     # Comprehensive QC report
│
├── qc/                         # Quality control metrics
│   ├── mosdepth/               # Coverage statistics per sample
│   ├── nanoplot/               # Read quality metrics (ONT)
│   └── cramino/                # CRAM-based QC (optional)
│
├── structural_variants/        # 🎯 PRIMARY OUTPUT: SV calls
│   ├── sniffles/               # Per-sample Sniffles VCFs
│   │   └── {sample}.sniffles.vcf.gz
│   ├── cutesv/                 # Per-sample CuteSV VCFs
│   │   └── {sample}.cutesv.vcf.gz
│   ├── svim/                   # Per-sample SVIM VCFs (if enabled)
│   │   └── {sample}.svim.vcf.gz
│   ├── merged/                 # Multi-caller consensus SVs
│   │   ├── {sample}.jasmine.vcf.gz
│   │   └── {sample}.survivor.vcf.gz
│   ├── annotated/              # AnnotSV annotations
│   │   └── {sample}.annotated.tsv
│   └── svanna/                 # Phenotype-prioritized SVs
│       └── {sample}.svanna.html
│
├── phasing/                    # Haplotype-resolved results
│   ├── haplotagged_bams/       # Phase-tagged alignments
│   │   └── {sample}.haplotagged.bam
│   ├── whatshap/               # Phasing statistics
│   │   └── {sample}.phased.vcf.gz
│   └── longphase/              # Alternative phasing
│       └── {sample}.longphase.vcf.gz
│
├── snv_calls/                  # SNVs (if --snv enabled)
│   ├── clair3/
│   │   └── {sample}.clair3.vcf.gz
│   └── deepvariant/
│       └── {sample}.deepvariant.vcf.gz
│
├── cnv_calls/                  # CNVs (if --cnv enabled)
│   └── spectre/
│       └── {sample}.cnv.vcf.gz
│
├── str_calls/                  # STRs (if --str enabled)
│   └── straglr/
│       └── {sample}.straglr.tsv
│
└── methylation/                # Methylation (if --methyl enabled, ONT only)
    └── modkit/
        └── {sample}.bedmethyl.gz
```

**Key output files:**

- **Merged SVs**: `structural_variants/merged/{sample}.jasmine.vcf.gz` (high-confidence consensus)
- **Annotated SVs**: `structural_variants/annotated/{sample}.annotated.tsv` (clinical interpretation)
- **QC Report**: `pipeline_info/multiqc_report.html` (overall quality assessment)
- **Phenotype-prioritized**: `structural_variants/svanna/{sample}.svanna.html` (ranked by phenotype match)

---

## Configuration Profiles

**Available Profiles:**

- test: Minimal test dataset
- docker: Use Docker containers
- singularity: Use Singularity containers

**Custom Configuration**

```bash
// custom.config
params {
    max_cpus = 16
    max_memory = '64.GB'
    outdir = '/scratch/results'
}

process {
    withName: 'CLAIR3' {
        cpus = 8
        memory = '32.GB'
    }
}
```

Run with:

```bash
nextflow run main.nf -c custom.config -profile docker
```

---

## Family-Based Analysis

### Trio/Family Configuration

For family-based SV analysis, provide pedigree information in your samplesheet:

```csv
sample,bam,bai,family,paternal_id,maternal_id,sex,phenotype,hpo_terms
child_001,child.bam,child.bam.bai,FAM001,father_001,mother_001,2,affected,"HP:0001250,HP:0002066,HP:0001263"
father_001,father.bam,father.bam.bai,FAM001,0,0,1,unaffected,
mother_001,mother.bam,mother.bam.bai,FAM001,0,0,2,unaffected,
```

**Sex encoding:** `1` = male, `2` = female, `0` = unknown
**Parental IDs:** Use `0` for founders (individuals with no parents in the study)

### De Novo SV Detection

Enable trio analysis to identify _de novo_ structural variants:

```bash
nextflow run nf-core/longraredisease \
    --input trio_samplesheet.csv \
    --trio_analysis true \
    --outdir trio_results \
    --fasta GRCh38.fasta \
    --sequencing_platform ont \
    -profile docker
```

The pipeline will:

1. ✅ Call SVs in each family member independently
2. ✅ Merge calls using JASMINE with family-aware parameters
3. ✅ Identify variants present in child but absent in parents
4. ✅ Filter based on read support and quality metrics

### Phenotype-Driven Prioritization (SVANNA)

When HPO terms are provided, SVANNA ranks SVs by phenotype relevance:

```bash
nextflow run nf-core/longraredisease \
    --input trio_samplesheet.csv \
    --trio_analysis true \
    --run_svanna true \
    --svanna_db /path/to/svanna/2302 \
    --outdir prioritized_results \
    --fasta GRCh38.fasta \
    --sequencing_platform ont \
    -profile docker
```

**Required:** Download SVANNA database from [Monarch Initiative](https://github.com/TheJacksonLaboratory/SvAnna)

**Output:** HTML report ranking SVs by:

- Overlap with disease-associated genes
- Regulatory impact predictions
- Phenotype similarity scores
- De novo status (if trio data available)

### Annotation with AnnotSV

Enable comprehensive SV annotation:

```bash
nextflow run nf-core/longraredisease \
    --input samplesheet.csv \
    --annotsv_db /path/to/AnnotSV_db \
    --outdir annotated_results \
    --fasta GRCh38.fasta \
    --sequencing_platform ont \
    -profile docker
```

**AnnotSV provides:**

- Gene overlap and functional impact
- ClinGen/ClinVar annotations
- DGV/gnomAD population frequencies
- Pathogenicity predictions (ACMG criteria)
- Regulatory element disruption

---

## Troubleshooting

### Common Issues

#### 1. Low SV Detection Rate

**Symptoms:** Fewer SVs than expected

**Solutions:**

```bash
# Lower read support threshold
--min_read_support 2

# Reduce minimum SV size
--min_sv_size 20

# Enable SVIM for better breakend detection
--run_svim true

# Lower consensus requirement
--jasmine_min_support 1
```

#### 2. High False Positive Rate

**Symptoms:** Many low-quality SV calls

**Solutions:**

```bash
# Increase read support
--min_read_support 5

# Require multiple caller agreement
--jasmine_min_support 2

# Increase minimum SV size
--min_sv_size 50
```

#### 3. Memory Issues

**Symptoms:** Process killed due to OOM

**Solutions:**

```bash
# Increase max memory
--max_memory 128.GB

# Reduce parallel processes
--max_cpus 16

# Use chromosome-based parallelization (automatic)
```

#### 4. Missing De Novo Variants

**Symptoms:** Expected _de novo_ variants not detected

**Checklist:**

- ✅ Ensure `--trio_analysis true` is set
- ✅ Verify pedigree information in samplesheet
- ✅ Check read coverage in all samples (≥30×)
- ✅ Review `structural_variants/merged/` for family calls
- ✅ Lower `--jasmine_min_support` if needed

#### 5. SVANNA Database Issues

**Symptoms:** SVANNA fails or produces no rankings

**Solutions:**

```bash
# Verify database path and version
ls -lh /path/to/svanna/2302

# Ensure HPO terms are valid (HP:XXXXXXX format)
# Check samplesheet for proper HPO term formatting

# Download latest SVANNA database:
wget https://storage.googleapis.com/svanna-db/svanna-data-2302.tar.gz
tar -xzf svanna-data-2302.tar.gz
```

### Performance Optimization

**For large cohorts (>10 samples):**

```bash
# Enable resource-efficient mode
--max_cpus 64
--max_memory 256.GB

# Use Singularity for better resource isolation
-profile singularity

# Enable work directory cleanup
-resume -with-dag flowchart.html
```

**For whole genome sequencing:**

- Expect 8-24 hours runtime (depending on coverage)
- Allocate 64-128GB RAM per sample for SV calling
- Use SSD storage for work directory (I/O intensive)

---

## Test Data

The pipeline includes test data for validation:

- Location: assets/test_data/
- Genome: Chromosome 22 subset
- Samples: Simulated nanopore data
- Runtime: ~10-15 minutes

---

## Getting Help

**Debugging Failed Runs:**

```bash
# Check Nextflow log for detailed errors
less .nextflow.log

# Resume from last successful step
nextflow run nf-core/longraredisease -resume

# Enable debug mode for verbose output
nextflow run nf-core/longraredisease --debug -profile docker
```

**Reporting Issues:**

When reporting issues, please include:

- Nextflow version (`nextflow -version`)
- Command used to run the pipeline
- Relevant error messages from `.nextflow.log`
- Sample metadata (anonymized if sensitive)
- System specifications (CPU, RAM, storage)

---

## Citation

If you use **nf-core/longraredisease** in your research, please cite:

> **nf-core/longraredisease: A Nextflow pipeline for long-read sequencing analysis in rare disease research** > _Citation to be added upon publication_

Additionally, please cite the tools used in your analysis:

**Core SV Tools:**

- **Sniffles2:** Sedlazeck et al. (2018) _Nature Methods_
- **CuteSV:** Jiang et al. (2020) _Genome Biology_
- **JASMINE:** Kirsche et al. (2023) _Nature Methods_
- **LongPhase:** Luo et al. (2023) _Nature Communications_
- **AnnotSV:** Geoffroy et al. (2018) _Bioinformatics_

**Optional Analysis Tools:**

- **SVANNA:** Danis et al. (2022) _AJHG_
- **Clair3:** Zheng et al. (2022) _Nature Computational Science_
- **Spectre:** Suvakov et al. (2021) _Genome Research_
- **Straglr:** Chin et al. (2023) _Genome Research_

---

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Make your changes following [nf-core guidelines](https://nf-co.re/developers/guidelines)
4. Test with `nextflow run . -profile test,docker`
5. Commit your changes (`git commit -m 'Add AmazingFeature'`)
6. Push to the branch (`git push origin feature/AmazingFeature`)
7. Open a Pull Request

**Please ensure:**

- ✅ Code follows nf-core style guidelines
- ✅ All tests pass successfully
- ✅ Documentation is updated accordingly
- ✅ Commit messages are descriptive

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

This pipeline was developed with support from [institution/funding sources]. We thank the nf-core community for infrastructure and best practices, and all tool developers whose software makes this pipeline possible.

---

**Pipeline Version:** 1.0.0
**Nextflow Version:** ≥25.10.4
**Last Updated:** 2024
