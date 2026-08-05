# 515Y–926R eASV Snakemake pipeline

This workflow is in progress and is not yet ready for production use.

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

The cloned repository intentionally does not include a `config/` folder. Create
and download the study-specific setup package from the pipeline tutorial, then
place its complete `config/` folder in the repository before running Snakemake.

## Run report

The final workflow target creates a self-contained report at:

```text
results/07-report/<studyName>.pipeline-report.html
```

The report summarizes the configuration, read retention, sample-level quality
control, domain composition, abundant taxa, and—when enabled—the
internal-standard results and figures.

For convenient downloading, the workflow also copies the report and formatted
data tables into the top-level `Results-Export/` folder. When internal standards
are enabled, the corrected ASV table is included there as well.
