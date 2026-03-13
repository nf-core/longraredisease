# nf-core/longraredisease: Output

## Introduction

This document describes the output produced by the pipeline, with a focus on **structural variant (SV) detection and analysis**. The pipeline is optimized for long-read sequencing data from Oxford Nanopore and PacBio platforms, providing comprehensive SV calling with multi-caller consensus and family-based analysis.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory and organized by sample ID.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and specializes in structural variant detection from long-read sequencing data. The core SV analysis workflow includes:

- [Quality Control](#quality-control) - Read QC and alignment statistics
- [Alignment](#alignment) - Read alignment to reference genome
- [Structural Variants](#structural-variants) - **Multi-caller SV detection, merging, and annotation** (PRIMARY FOCUS)
- [Phasing](#phasing) - Haplotype phasing for improved SV calling
- [Family Analysis](#family-analysis) - Trio-based SV calling and inheritance patterns
- [SV Prioritization](#sv-prioritization) - Phenotype-driven variant ranking
- [MultiQC](#multiqc) - Aggregate report with QC metrics
- [Pipeline information](#pipeline-information) - Execution reports and metadata

**Optional analyses** (disabled by default in SV-focused mode):

- Single Nucleotide Variants - SNV/indel calling
- Copy Number Variants - CNV detection
- Short Tandem Repeats - STR genotyping
- Methylation - DNA methylation calling

## SV Detection Strategy

The pipeline uses a **multi-caller consensus approach** for high-confidence SV detection:

```
┌─────────────────────────────────────────────────────────────┐
│                      ALIGNED READS                          │
│                  (minimap2 + samtools)                      │
└────────────┬────────────────────────────────────────────────┘
             │
             ├──────────────────────────────────────────┐
             │                                          │
             ▼                                          ▼
    ┌──────────────────┐                      ┌─────────────────┐
    │  PHASING         │                      │  SNV CALLING    │
    │  (LongPhase)     │                      │  (Clair3 / DV)  │
    └────────┬─────────┘                      └─────────────────┘
             │
             ▼
    ┌──────────────────┐
    │  HAPLOTAGGING    │
    │  (BAM + HP tags) │
    └────────┬─────────┘
             │
             ▼
    ┌──────────────────────────────────────────────────────────┐
    │              MULTI-CALLER SV DETECTION                   │
    │                                                           │
    │  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌────────────┐ │
    │  │Sniffles │  │ CuteSV  │  │  SVIM   │  │ Trio Mode  │ │
    │  │(primary)│  │(option) │  │(option) │  │ (if fam)   │ │
    │  └────┬────┘  └────┬────┘  └────┬────┘  └──────┬─────┘ │
    └───────┼────────────┼────────────┼───────────────┼───────┘
            │            │            │               │
            └────────────┴────────────┴───────────────┘
                         │
                         ▼
            ┌────────────────────────┐
            │  MERGING & CONSENSUS   │
            │      (JASMINE)         │
            └───────────┬────────────┘
                        │
                        ▼
            ┌────────────────────────┐
            │  ANNOTATION            │
            │  (AnnotSV)             │
            └───────────┬────────────┘
                        │
                        ▼
            ┌────────────────────────┐
            │  PRIORITIZATION        │
            │  (SVANNA + HPO terms)  │
            └────────────────────────┘
```

**Key features:**

- **Sniffles** (primary caller) - Fast, accurate, phase-aware SV detection
- **CuteSV** (optional) - Complementary caller for DEL/INS/DUP
- **SVIM** (optional, `--run_svim true`) - Recommended for complex SVs and breakends (BND)
- **JASMINE merging** - Consensus calling across all enabled callers
- **Trio-based calling** - When `--trio_analysis true`, uses family information for de novo detection
- **AnnotSV** - Clinical annotation with disease databases
- **SVANNA** - HPO-based phenotype matching for clinical prioritization

## Output directory structure

```
results/
├── <sample_id>/                    # Per-sample results directory
│   ├── pedigree_file/             # Family pedigree information (if trio analysis)
│   ├── mapped_bam/                # Aligned BAM files
│   ├── nanoplot_qc/               # NanoPlot QC reports
│   ├── bam_stats/                 # BAM statistics (flagstat, stats, idxstats)
│   ├── mosdepth/                  # Coverage analysis
│   │
│   ├── sniffles/                  # Sniffles SV calls (primary SV caller)
│   ├── cutesv/                    # CuteSV SV calls (optional)
│   ├── svim/                      # SVIM SV calls (optional, better for BNDs)
│   ├── merged_sv/                 # Merged multi-caller SV consensus (JASMINE)
│   ├── filtered_pass_sv/          # PASS-filtered SV calls
│   ├── downsampled_sv/            # Downsampled SV calls (if coverage-based filtering)
│   ├── annotsv_sniffles/          # AnnotSV annotations for Sniffles calls
│   ├── annotsv_svim/              # AnnotSV annotations for SVIM calls (if run)
│   ├── svanna/                    # Phenotype-driven SV prioritization (if HPO terms provided)
│   │
│   ├── sniffles_trio/             # Trio-based joint SV calling (if family analysis)
│   ├── unphased_sniffles/         # Unphased Sniffles calls (if haplotagging disabled)
│   │
│   ├── longphase/                 # Phased variants for SV calling
│   ├── haplotagged_bam/           # Haplotagged BAM for phased SV detection
│   │
│   └── unify_vcf/                 # Unified VCF output (optional)
│
├── multiqc/                        # Aggregate MultiQC report
├── pipeline_info/                  # Pipeline execution metadata
├── annotsv_db/                     # AnnotSV annotation database (if downloaded)
└── test_family/                    # Family-level outputs (if trio analysis enabled)
    └── sniffles_trio/             # Joint-called family SV variants
```

## Quality Control

### NanoPlot

<details markdown="1">
<summary>Output files</summary>

- `qc/nanoplot/<sample_id>/`
  - `NanoPlot-report.html`: Interactive HTML report with read statistics and plots
  - `NanoStats.txt`: Summary statistics in text format
  - `*_HistogramReadlength.png`: Read length distribution
  - `*_LengthvsQualityScatterPlot.png`: Read length vs quality scatter plot
  - `*_Weighted_HistogramReadlength.png`: Weighted read length histogram
  - `*_Weighted_LogTransformed_HistogramReadlength.png`: Log-transformed weighted histogram

</details>

[NanoPlot](https://github.com/wdecoster/NanoPlot) generates comprehensive quality control plots and statistics for long-read sequencing data. The reports include:

- Read length distributions
- Read quality distributions
- Sequencing throughput over time
- Read length vs quality scatter plots

### Mosdepth

<details markdown="1">
<summary>Output files</summary>

- `qc/mosdepth/<sample_id>/`
  - `*.mosdepth.global.dist.txt`: Global coverage distribution
  - `*.mosdepth.region.dist.txt`: Per-region coverage distribution (if BED file provided)
  - `*.mosdepth.summary.txt`: Summary statistics for coverage
  - `*.per-base.bed.gz`: Per-base depth (if enabled)
  - `*.regions.bed.gz`: Per-region depth (if BED file provided)

</details>

[Mosdepth](https://github.com/brentp/mosdepth) provides fast BAM/CRAM depth calculation with coverage statistics at different levels. Key metrics include:

- Mean coverage across the genome
- Coverage uniformity
- Percentage of bases covered at different thresholds (e.g., 10x, 20x, 30x)

## Alignment

### Minimap2 / Winnowmap

<details markdown="1">
<summary>Output files</summary>

- `alignment/<sample_id>/`
  - `*.sorted.bam`: Coordinate-sorted aligned BAM file
  - `*.sorted.bam.bai`: BAM index file
  - `*.flagstat`: Samtools flagstat alignment statistics

</details>

The pipeline aligns reads using either [Minimap2](https://github.com/lh3/minimap2) (default) or [Winnowmap](https://github.com/marbl/Winnowmap) for repetitive regions. Alignment parameters are automatically optimized based on the sequencing platform (ONT, PacBio HiFi, PacBio CLR).

## Structural Variants

The pipeline implements a **multi-caller SV detection strategy** to maximize sensitivity and specificity. SVs are detected independently by up to three complementary tools, then merged to create a high-confidence consensus callset.

### Sniffles (Primary SV Caller)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/sniffles/`
  - `<sample_id>_sniffles.vcf.gz`: Structural variant calls in VCF format
  - `<sample_id>_sniffles.vcf.gz.tbi`: VCF index file
  - `<sample_id>_sniffles.snf`: Sniffles internal file for joint calling (if trio analysis enabled)

</details>

[Sniffles2](https://github.com/fritzsedlazeck/Sniffles) is the **primary SV caller**, detecting structural variants from long-read alignments with high sensitivity and precision.

**Configuration:**

- Minimum mapping quality: `--sniffles_min_mapq` (default: 10)
- Phasing support: Enabled when `--haplotag_bam` is true
- Tandem repeat masking: Uses `sniffles_tandem_file` to improve specificity

**SV types detected:**

- **DEL** - Deletions
- **INS** - Insertions
- **DUP** - Duplications
- **INV** - Inversions
- **BND** - Translocations/breakends

### CuteSV (Complementary Caller)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/cutesv/`
  - `<sample_id>_cutesv_final.vcf.gz`: Structural variant calls with improved read support
  - `<sample_id>_cutesv_final.vcf.gz.tbi`: VCF index

</details>

[CuteSV](https://github.com/tjiangHIT/cuteSV) provides complementary SV calling with optimized detection of **insertions and complex variants**. It uses a different algorithmic approach than Sniffles, improving overall sensitivity when calls are merged.

**Configuration:**

- Minimum mapping quality: `--cutesv_min_mapq` (default: 10)
- Genotyping: Enabled by default
- Read support adjustment: Processed by `RE2SUPPORT` module to standardize support metrics

**Enabled when:** `--merge_sv true`

### SVIM (Breakend Specialist)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/svim/`
  - `<sample_id>_svim_final.vcf.gz`: Structural variant calls
  - `<sample_id>_svim_final.vcf.gz.tbi`: VCF index

</details>

[SVIM](https://github.com/eldariont/svim) uses alignment signatures to detect SVs with particular strength in identifying **complex breakends (BND)** and nested variants.

**Configuration:**

- Minimum mapping quality: `--svim_min_mapq` (default: 10)
- Interspersed duplications reported as insertions
- Tandem duplications reported as insertions

**Enabled when:** `--run_svim true` (recommended for comprehensive breakend detection)

### Merged Multi-Caller SV Consensus (JASMINE)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/merged_sv/`
  - `<sample_id>.final.vcf.gz`: Merged SV calls from multiple callers
  - `<sample_id>.final.vcf.gz.tbi`: VCF index

</details>

[JASMINE](https://github.com/mkirsche/Jasmine) merges SV calls from multiple callers (Sniffles, CuteSV, SVIM) to create a **high-confidence consensus callset**. This approach:

✅ **Reduces false positives** - SVs supported by multiple callers are more likely to be true variants
✅ **Improves genotyping** - Combines evidence across callers for better genotype accuracy
✅ **Retains caller-specific calls** - High-quality calls from individual callers are preserved

**Configuration:**

- Minimum caller support: `--min_caller_support` (default: 2) - Requires at least N callers to agree
- Variant ID tracking: `--keep_var_ids` preserves original caller IDs
- Genotype output: Enabled by default

**Enabled when:** `--merge_sv true`

### Filtered SV Calls

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/filtered_pass_sv/`
  - `<sample_id>.<caller>.filtered.vcf.gz`: PASS-only filtered variants
  - `<sample_id>.<caller>.filtered.vcf.gz.tbi`: VCF index

</details>

When `--filter_pass_sv true`, only variants with `FILTER=PASS` are retained. This removes low-quality calls flagged by the SV callers.

### Downsampled SV Calls

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/downsampled_sv/`
  - `<sample_id>.<caller>.downsampled.vcf.gz`: Coverage-based filtered variants
  - `<sample_id>.<caller>.downsampled.vcf.gz.tbi`: VCF index

</details>

When `--downsample_sv true` and `--coverage_bed` is provided, SVs are filtered based on coverage thresholds to remove calls from low-coverage or repetitive regions.

**Configuration:**

- Minimum read support: `--min_read_support` (default: 'auto')
- Minimum read support limit: `--min_read_support_limit` (default: 2)

### Annotated SVs (AnnotSV)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotsv_sniffles/`

  - `<sample_id>_annotsv.tsv`: Comprehensive tab-separated annotation file
  - `<sample_id>_annotsv.unannotated.tsv`: Unannotated variant subset

- `<sample_id>/annotsv_svim/` (if `--run_svim true`)
  - `<sample_id>_annotsv.tsv`: SVIM-specific annotations

</details>

[AnnotSV](https://github.com/lgmgeo/AnnotSV) provides **comprehensive SV annotation** with clinical and population genetics databases.

**Annotations include:**

- 🧬 **Gene overlap** - Genes affected by SV breakpoints
- 📊 **Population frequencies** - gnomAD-SV, DGV, DDD databases
- ⚕️ **Clinical significance** - ClinVar pathogenic variants, OMIM diseases
- 🎯 **Functional predictions** - Haploinsufficiency scores, regulatory regions
- 👨‍👩‍👧 **Pathogenicity ranking** - ACMG-based classification (`--rankfiltering`)

**Configuration:**

- Genome build: GRCh38 (hardcoded)
- Rank filtering: `--rankfiltering` (default: '3-5,NA') - Reports pathogenic to VUS variants
- HPO terms: Automatically used if provided in samplesheet for phenotype-driven ranking
- Custom annotations: `--annotsv_annotations` path to additional annotation sources

**Enabled when:** `--annotate_sv true`

## Family Analysis

### Pedigree File

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/pedigree_file/`
  - `<family_id>.ped`: PED format pedigree file for trio analysis

</details>

When `--trio_analysis true`, the pipeline generates a pedigree file from the samplesheet family relationships. This file is used for:

- Joint SV calling with family information
- Mendelian inheritance analysis
- De novo variant detection

### Trio-based Joint SV Calling

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/sniffles_trio/`

  - `<sample_id>_trio.vcf.gz`: Joint-called family SV variants
  - `<sample_id>_trio.vcf.gz.tbi`: VCF index

- `test_family/sniffles_trio/`
  - Family-level joint-called variants (organized by family ID)

</details>

When `--trio_analysis true`, Sniffles performs **joint SV calling** across all family members simultaneously. This approach:

✅ **Improves accuracy** - Uses family information to refine genotype calls
✅ **Detects de novo variants** - Identifies SVs present in child but not in parents
✅ **Validates inheritance** - Confirms Mendelian inheritance patterns
✅ **Increases sensitivity** - Detects low-coverage SVs supported by family structure

**Requirements:**

- Samplesheet must include `family`, `paternal_id`, and `maternal_id` columns
- All family members must be processed in the same run

### Mendelian Inheritance Analysis (RTG Tools)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/rtg_mendelian_sv/`

  - `<sample_id>_mendelian.vcf.gz`: Variants following Mendelian inheritance
  - `<sample_id>_mendelian_stats.txt`: Mendelian consistency statistics

- `<sample_id>/rtg_violations_sv/`
  - `<sample_id>_violations.vcf.gz`: Mendelian violations (potential de novo or errors)
  - `<sample_id>_violations_stats.txt`: Violation statistics

</details>

RTG Tools analyzes family-based SV calls to identify:

- **Mendelian-consistent variants** - SVs following expected inheritance patterns
- **Mendelian violations** - Candidate de novo SVs or genotyping errors
- **Inheritance patterns** - Dominant, recessive, compound heterozygous

## SV Prioritization

### SVANNA (Phenotype-driven Prioritization)

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/svanna/`
  - `<sample_id>_svanna.html`: Interactive HTML report with prioritized variants
  - `<sample_id>_svanna.tsv`: Tab-separated prioritization results
  - `<sample_id>_svanna.vcf`: Annotated VCF with SVANNA scores

</details>

[SVANNA](https://github.com/TheJacksonLaboratory/SvAnna) performs **phenotype-driven SV prioritization** using Human Phenotype Ontology (HPO) terms. It ranks variants based on:

🎯 **Phenotype relevance** - Matches SV-affected genes to patient phenotypes
🧬 **Functional impact** - Predicts effect on gene function
📚 **Disease associations** - Links to known disease-gene relationships
🔬 **Variant characteristics** - Considers SV type, size, and breakpoint precision

**Configuration:**

- HPO terms: Provided in samplesheet `hpo_terms` column (e.g., "HP:0001250,HP:0002066")
- Database: `--svanna_db` path to SVANNA database
- Output format: `--output_format` (html, tsv, vcf)

**Enabled when:** `--run_svanna true` AND HPO terms are provided

**Use case:** Essential for clinical interpretation of rare disease cases where patient phenotype can guide variant prioritization.

## Single Nucleotide Variants

_These analyses are **optional** and disabled by default in SV-focused mode. Enable with `--snv true`._

### Clair3

<details markdown="1">
<summary>Output files</summary>

- `snv/clair3/<sample_id>/`
  - `*.clair3.vcf.gz`: SNV and indel calls
  - `*.clair3.vcf.gz.tbi`: VCF index
  - `*_pileup.vcf.gz`: Pileup variant calls
  - `*_full_alignment.vcf.gz`: Full alignment variant calls

</details>

[Clair3](https://github.com/HKU-BAL/Clair3) uses deep learning for accurate SNV and small indel calling from long reads. It provides:

- High-quality SNV calls with genotype likelihoods
- Small indel detection
- Support for different basecalling models

### DeepVariant

<details markdown="1">
<summary>Output files</summary>

- `snv/deepvariant/<sample_id>/`
  - `*.deepvariant.vcf.gz`: SNV and indel calls
  - `*.deepvariant.vcf.gz.tbi`: VCF index
  - `*.deepvariant.g.vcf.gz`: gVCF file for joint calling

</details>

[DeepVariant](https://github.com/google/deepvariant) provides complementary SNV calling using convolutional neural networks. It can be used alongside Clair3 for higher confidence calls.

### Joint-genotyped SNVs (Trio analysis)

<details markdown="1">
<summary>Output files</summary>

- `snv/joint_genotyped/<family_id>/`
  - `*.joint.vcf.gz`: Joint-genotyped VCF for all family members
  - `*.joint.vcf.gz.tbi`: VCF index
  - `*.denovo.vcf.gz`: Predicted de novo variants (if trio)
  - `*.inheritance.txt`: Inheritance pattern analysis

</details>

When trio analysis is enabled, SNV calls from family members are joint-genotyped using [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools) to improve accuracy and identify de novo mutations.

### Annotated SNVs

<details markdown="1">
<summary>Output files</summary>

- `snv/annotated/<sample_id>/`
  - `*.snpeff.vcf.gz`: SNV calls annotated with SnpEff
  - `*.snpeff_summary.html`: SnpEff annotation summary
  - `*.vep.vcf.gz`: SNV calls annotated with VEP (if enabled)

</details>

SNVs are functionally annotated to predict their impact on genes and proteins, including:

- Variant consequences (missense, nonsense, splice site, etc.)
- Gene and transcript information
- Conservation scores
- Population frequencies
- Pathogenicity predictions (SIFT, PolyPhen, etc.)

## Copy Number Variants

_Optional analysis - disabled by default in SV-focused mode. Enable with `--cnv true`._

### Spectre

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/spectre/`
  - `<sample_id>.wf_cnv.vcf.gz`: CNV calls in VCF format
  - `<sample_id>.wf_cnv.bed`: CNV calls in BED format

</details>

[Spectre](https://github.com/fritzsedlazeck/Spectre) detects copy number variants from long-read sequencing data by analyzing coverage and structural variant signatures. Spectre is recommended for ONT data.

### HiFiCNV

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/hificnv/`
  - `<sample_id>.hificnv.vcf.gz`: CNV calls
  - `<sample_id>.hificnv.bed`: CNV regions

</details>

[HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) provides CNV detection optimized for PacBio HiFi data. Enable with `--cnv_hificnv true`.

## Short Tandem Repeats

_Optional analysis - disabled by default in SV-focused mode. Enable with `--str true`._

### Straglr

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/straglr/`

  - `<sample_id>_straglr_sorted.vcf.gz`: STR genotypes in VCF format
  - `<sample_id>_straglr.tsv`: STR genotypes in TSV format

- `<sample_id>/stranger/`
  - Annotated STR calls with pathogenic expansion information

</details>

[Straglr](https://github.com/bcgsc/straglr) genotypes short tandem repeats and detects expansions associated with genetic diseases such as Huntington's disease, fragile X syndrome, and myotonic dystrophy.

## Methylation

_Optional analysis - disabled by default in SV-focused mode. Enable with `--methyl true`. **ONT-specific** - requires base modification data._

### Modkit

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/methyl/`

  - `<sample_id>_methyl.bed`: Methylation calls in BED format
  - `<sample_id>_cpg.methyl.bed`: CpG-specific methylation calls

- `<sample_id>/methyl_bedgraph/`
  - `*.bedgraph`: Methylation bedgraph for genome browser visualization

</details>

[Modkit](https://github.com/nanoporetech/modkit) extracts modified base calls from Oxford Nanopore sequencing data with base modification models. Detects 5mC, 6mA, and other base modifications.

## Phasing

Phasing is **critical for accurate SV detection** in long-read data. The pipeline uses LongPhase to assign variants to haplotypes, which enables:

- Phase-aware SV calling with Sniffles
- Haplotype-resolved SV visualization
- Improved SV genotyping accuracy
- Detection of compound heterozygous SVs

### LongPhase

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/longphase/`
  - `<sample_id>_longphase.vcf.gz`: Phased variant calls
  - `<sample_id>_longphase.vcf.gz.tbi`: VCF index
  - `<sample_id>_longphase_stats.txt`: Phase block statistics

</details>

[LongPhase](https://github.com/twolinin/LongPhase) performs haplotype phasing using long-read information. It provides:

- **Long phase blocks** - Leverages long-read spanning multiple variants
- **Phased SNVs** - Used as markers for haplotagging reads
- **Platform-specific optimization** - Automatic selection of `--ont` or `--pb` mode

**Configuration:**

- Platform: Automatically set based on `--sequencing_platform`
- Input variants: Uses Clair3 or DeepVariant SNV calls for phasing backbone

### Haplotagged BAMs

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/haplotagged_bam/`
  - `<sample_id>_haplotagged.bam`: BAM file with haplotype tags (HP tag)
  - `<sample_id>_haplotagged.bam.bai`: BAM index

</details>

When `--haplotag_bam true`, reads are tagged with their assigned haplotype (HP:1 or HP:2). **Haplotagged BAMs are used by Sniffles** for:

✅ **Phase-aware SV calling** - Sniffles uses HP tags to assign SVs to haplotypes
✅ **Improved specificity** - Reduces false positives from haplotype-switching artifacts
✅ **Compound heterozygote detection** - Identifies SVs on opposite haplotypes
✅ **Visualization** - Enables haplotype-colored IGV tracks

**Recommended:** Enable haplotagging for all SV-focused analyses

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Standalone HTML report viewable in web browser
  - `multiqc_data/`: Directory with parsed statistics from all tools
  - `multiqc_plots/`: Static images from the report

</details>

[MultiQC](http://multiqc.info) generates a single comprehensive HTML report summarizing quality control metrics and analysis results from all samples and tools in the pipeline. The report includes:

- Alignment statistics
- Coverage metrics
- Variant calling summaries
- Tool version information
- Run statistics

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report with resource usage
  - `execution_timeline.html`: Timeline of process execution
  - `execution_trace.txt`: Detailed trace of all processes
  - `pipeline_dag.svg`: Pipeline directed acyclic graph visualization
  - `software_versions.yml`: Software versions used in the pipeline
  - `params.json`: Parameters used for the pipeline run

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides detailed reports about pipeline execution for troubleshooting, resource optimization, and reproducibility.

---

## Quick Reference: Key Output Files

### For Clinical SV Interpretation

| Analysis Stage             | File Location                                           | Description                                  |
| -------------------------- | ------------------------------------------------------- | -------------------------------------------- |
| **Primary SV Calls**       | `<sample>/sniffles/<sample>_sniffles.vcf.gz`            | Phased SVs from Sniffles (primary caller)    |
| **Multi-caller Consensus** | `<sample>/merged_sv/<sample>_merged.vcf.gz`             | High-confidence SVs from JASMINE merging     |
| **Clinical Annotations**   | `<sample>/annotsv_sniffles/<sample>.annotated.tsv`      | AnnotSV annotations with disease information |
| **Prioritized SVs**        | `<sample>/svanna/<sample>_svanna.html`                  | Phenotype-ranked SVs (if HPO terms provided) |
| **Trio Analysis**          | `<sample>/sniffles_trio/<sample>_trio.vcf.gz`           | Family-based SV calls with de novo detection |
| **Mendelian Violations**   | `<sample>/rtg_violations_sv/<sample>_violations.vcf.gz` | Candidate de novo SVs                        |

### For Quality Assessment

| QC Metric            | File Location                                     | What to Check                                 |
| -------------------- | ------------------------------------------------- | --------------------------------------------- |
| **Read QC**          | `<sample>/nanoplot_qc/NanoPlot-report.html`       | Read quality, N50, coverage distribution      |
| **Alignment**        | `<sample>/bam_stats/<sample>.flagstat.txt`        | Mapping rate (should be >95%)                 |
| **Coverage**         | `<sample>/mosdepth/<sample>.mosdepth.summary.txt` | Mean coverage (recommend ≥30x for SV calling) |
| **Pipeline Summary** | `multiqc/multiqc_report.html`                     | Aggregate QC across all samples               |

### For Downstream Analysis

| Analysis Type                    | File to Use                                               | Format                                |
| -------------------------------- | --------------------------------------------------------- | ------------------------------------- |
| **Genome Browser Visualization** | `<sample>/haplotagged_bam/<sample>_haplotagged.bam`       | Haplotype-colored BAM                 |
| **Variant Filtering**            | `<sample>/filtered_pass_sv/<sample>_filtered_pass.vcf.gz` | VCF with PASS filter                  |
| **Database Upload**              | `<sample>/annotsv_sniffles/<sample>.annotated.tsv`        | Tab-delimited with all annotations    |
| **Further Annotation**           | `<sample>/merged_sv/<sample>_merged.vcf.gz`               | Clean merged VCF for custom pipelines |

### Common Analysis Patterns

#### Single Sample Discovery

```
1. Check QC: multiqc/multiqc_report.html
2. Review primary calls: <sample>/sniffles/<sample>_sniffles.vcf.gz
3. Explore annotations: <sample>/annotsv_sniffles/<sample>.annotated.tsv
4. Visualize candidates: Load <sample>/haplotagged_bam/<sample>_haplotagged.bam in IGV
```

#### Trio/Family Analysis

```
1. Verify pedigree: <sample>/pedigree_file/<family_id>.ped
2. Review joint calls: <sample>/sniffles_trio/<sample>_trio.vcf.gz
3. Find de novo variants: <sample>/rtg_violations_sv/<sample>_violations.vcf.gz
4. Prioritize with phenotype: <sample>/svanna/<sample>_svanna.html
```

#### High-confidence Rare Disease Diagnosis

```
1. Multi-caller consensus: <sample>/merged_sv/<sample>_merged.vcf.gz
2. Phenotype prioritization: <sample>/svanna/<sample>_svanna.html (if HPO terms)
3. Review top candidates in: <sample>/annotsv_sniffles/<sample>.annotated.tsv
4. Visual validation: IGV with haplotagged BAM
5. Check inheritance: RTG mendelian analysis (if trio)
```
