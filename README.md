<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-longraredisease_logo_dark.png">
    <img alt="nf-core/longraredisease" src="docs/images/nf-core-longraredisease_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/longraredisease)
[![GitHub Actions CI Status](https://github.com/nf-core/longraredisease/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/longraredisease/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/longraredisease/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/longraredisease/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/longraredisease/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.20935122-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.20935122)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/longraredisease)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23longraredisease-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/longraredisease)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/longraredisease** is a specialized bioinformatics pipeline for **structural variant (SV) detection and clinical interpretation** from long-read sequencing data (Oxford Nanopore and PacBio). Designed for rare disease diagnostics, it delivers high-confidence variant discovery through multi-caller consensus, family-based analysis, and phenotype-driven prioritization.

![Long-read sequencing pipeline](docs/images/longraredisease_pipeline.png)

The pipeline supports:

- **Multi-caller SV consensus** — Sniffles, CuteSV, SVIM with JASMINE merging
- **Phase-aware calling** — Haplotype-resolved SV detection using LongPhase
- **Family analysis** — Trio-based joint calling and de novo variant detection
- **Clinical annotation** — AnnotSV with disease database integration
- **Phenotype prioritization** — SVANNA-based ranking using HPO terms
- **Optional analyses** — SNVs (Clair3/DeepVariant), CNVs (Spectre/HiFiCNV), STRs (Straglr), Methylation (Modkit)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data:

```csv title="samplesheet.csv"
sample,file_path,hpo_terms,sex,phenotype,family_id,maternal_id,paternal_id
sample1,/path/to/sample1.bam,HP:0002721;HP:0002110,1,2,,,
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/longraredisease \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --fasta reference.fasta \
    --sequencing_platform ont
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided with the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details on pipeline usage and parameters, see [docs/usage.md](docs/usage.md).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/longraredisease/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to [docs/output.md](docs/output.md).

## Credits

I thank the following people for their contributions and guidance to the development of the pipeline:

The [nf-core](https://nf-co.re) team, and especially Friederike Hanssen, Ken Brewer, Nicolas Vannieuwkerke and Maxime U Garcia for their support and guidance in developing this pipeline.

I also thank the clinical scientists Chipo Mashayamombe-Wolfgarten, Hannah Titheradge, and Lorraine Hartles-Spencer for their invaluable clinical input and expertise. I would also like to thank Professor Andrew Beggs for his clinical guidance and expertise.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#longraredisease` channel](https://nfcore.slack.com/channels/longraredisease) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/longraredisease for your analysis, please cite it using the following doi: [10.5281/zenodo.20935122](https://doi.org/10.5281/zenodo.20935122)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
