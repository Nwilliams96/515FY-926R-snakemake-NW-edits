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

## USC CARC runner

`run_snakemake_USC_CARC_only.sh` requests one eight-core, 128 GB CARC job and
passes the Slurm CPU count to Snakemake. BBSplit, DADA2, and taxonomic
classification use up to eight threads. Snakemake also receives a 120 GB
memory budget, while the memory-heavy BBSplit and DADA2 rules declare realistic
per-job requirements so they are not run concurrently on the same allocation.
This USC-specific script is optional; other systems can use `run_snakemake.sh`
or provide allocation settings appropriate for their scheduler.

From the project directory, submit the CARC job with:

```bash
bash run_snakemake_USC_CARC_only.sh
```

The script reads `studyName` from `config/config.yml`, submits itself through
Slurm, and labels the queue entry and log files as `<studyName>`. Calling
the file with `sbatch` still works, but the log filename initially uses the
generic `pipeline` label because Slurm opens it before the script runs.

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
