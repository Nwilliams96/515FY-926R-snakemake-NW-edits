# 515Y–926R eASV Snakemake pipeline

This workflow is in progress and is not yet ready for production use.

Use the [interactive pipeline tutorial and configuration builder](https://www.nathanlrwilliams.com/eASV-Pipeline-Tutorial/)
to prepare a study-specific setup package. The tutorial's HTTPS clone command
downloads the current version of this repository's `main` branch.

## Requirements

Install [Git](https://git-scm.com/downloads) and
[Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/),
then create an isolated Snakemake controller environment following the
[official Snakemake installation guidance](https://snakemake.readthedocs.io/en/latest/getting_started/installation.html):

```bash
conda create --name snakemake --channel conda-forge --channel bioconda snakemake
conda activate snakemake
snakemake --version
```

Run this workflow with `--use-conda`. On the first run, Snakemake downloads and
creates the rule-specific QIIME 2, R, BLAST, BBMap, and utility environments.
Subsequent runs reuse them.

When internal-standard correction is enabled, `intstds` is an ordered list and
may contain one or more unique standard IDs:

```yaml
use_internal_standards: true
intstds:
  - "ISD_1"
  - "ISD_2"
```

Each ID requires a matching `<ID>_ng` column in `config/samples.tsv` and a
matching row in `config/internal_stds.tsv`. Older configs that use the
`intstd1`/`intstd2`/`intstd3` mapping remain supported.

The cloned repository intentionally does not include a `config/` folder. Create
and download the study-specific setup package from the pipeline tutorial, then
place its complete `config/` folder in the repository before running Snakemake.

## DADA2 controls

The tutorial exposes separate DADA2 settings for the paired 16S path and the
concatenated, single-end 18S path. Generated configs contain a `dada2` block
covering expected-error filtering, quality truncation, paired-read overlap,
pooling, chimera detection, parent abundance, and error-model training reads.
The displayed presets reproduce the pipeline's historical settings. Configs
created before this block was introduced remain supported through matching
workflow defaults.

## Run report

The final workflow target creates a self-contained report at:

```text
results/07-report/<studyName>.pipeline-report.html
```

The report summarizes the configuration, read retention, sample-level quality
control, domain composition, 16S/18S/chloroplast/mitochondrial assignments,
unassigned sequences, abundant taxa, and—when enabled—the
internal-standard results and figures.

For convenient downloading, the workflow also copies the report and formatted
data tables into the top-level `Results-Export/` folder. When internal standards
are enabled, the corrected ASV table is included there as well.
