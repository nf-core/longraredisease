# nf-core/longraredisease: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0 - [2026-06-26]

Initial release of nf-core/longraredisease, created with the [nf-core](https://nf-co.re/) template.

This release provides an end-to-end pipeline for comprehensive variant detection — structural variants (SV), single-nucleotide variants (SNV), copy-number variants (CNV), short tandem repeats (STR) and methylation — and clinical interpretation from long-read sequencing data (Oxford Nanopore and PacBio), aimed at rare disease diagnostics.

### `Added`

- Read alignment to a reference genome with `minimap2` (or `winnowmap`), and BAM processing/QC with `samtools` (sort, index, merge, stats, flagstat, idxstats).
- Multi-caller structural variant discovery using Sniffles, CuteSV and SVIM, with consensus merging via JASMINE.
- Haplotype-aware (phased) SV calling using LongPhase (`phase` and `haplotag`).
- Family/trio-based joint calling and de novo SV detection.
- Clinical SV annotation with AnnotSV (including annotation database installation) and SnpEff.
- Phenotype-driven SV prioritisation with SVANNA using HPO terms.
- Optional single-nucleotide variant calling with Clair3 and DeepVariant.
- Optional copy-number variant calling (Spectre/HiFiCNV), short-tandem-repeat genotyping (TRGT/Straglr) and methylation pileups with Modkit.
- VCF normalisation, filtering, merging and indexing with `bcftools` and `tabix`.
- Per-run reporting aggregated with MultiQC.

### `Fixed`

### `Dependencies`

### `Deprecated`
