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

Run this workflow with `--use-conda`. The included launch scripts also pass a
shared `--conda-prefix`, taken from `conda_envs_dir` in `config/config.yml`.
Snakemake identifies rule environments from their environment definitions: it
reuses matching completed environments and creates only those that are missing
or whose definitions have changed. Keep this directory outside each analysis
clone so new projects can share the same QIIME 2, R, BLAST, BBMap, and utility
environments. For example:

```yaml
conda_envs_dir: "../eASV-conda-envs"
```

Set `SNAKEMAKE_CONDA_PREFIX` before running if a computer or HPC needs to
override that configured location. Older projects without `conda_envs_dir`
reuse their existing `.snakemake/conda` directory when present; otherwise the
scripts default to `../eASV-conda-envs`.

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

Internal-standard method tables use the configured IDs in their filenames. For
example, standards named `BP` and `DR` produce `asv_table_BP_recovery_ratio.tsv`
and `asv_table_mean_BP_and_DR_recovery_ratio.tsv`, so the correction represented
by every single-standard and combination output is explicit. The combined long
table similarly includes all configured IDs, for example
`<study>.BP_DR.ISD_corrected_long_data.tsv`.

The required starting-pool concentration file is named
`config/prok_and_euk_SSU_amplicon_concentrations.tsv`. Updated workflows still
accept an existing `config/bioanalyzer.tsv` when the new file is absent, allowing
older study folders to finish after a pipeline update.

Additional user-defined columns may be added to `config/samples.tsv`. They are
preserved in the generated QIIME metadata and are available as filters in the
HTML report; they do not alter read processing.

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
sbatch run_snakemake_USC_CARC_only.sh
```

The script reads `studyName` from `config/config.yml`, submits itself through
Slurm, and labels the queue entry and log files as `<studyName>`. Calling
the file with `sbatch` still works, but the log filename initially uses the
generic `pipeline` label because Slurm opens it before the script runs.
Before starting Snakemake, the CARC runner also detects config files whose
timestamps are ahead of the compute-node clock and resets only those timestamps.
This prevents clock-skew failures after moving ZIP packages between computers.
The runner gives every Slurm job a package cache under CARC's job-specific
`TMPDIR` while retaining completed rule environments in the configured shared
`conda_envs_dir`. This prevents simultaneous jobs on different nodes from
contending for Conda metadata in `~/.conda/pkgs` and avoids stale-file-handle
failures during environment creation.

When databases must be built, SILVA and PR2 preparation can run concurrently.
BBSplit indexing and both primer-extraction steps use up to eight cores. QIIME's
naive-Bayes classifier-training action does not expose a worker-count option, so
each classifier-training command remains single-core while the independent
SILVA and PR2 branches may still run in parallel.

## Run report

The final workflow target creates a self-contained report at:

```text
results/07-report/<studyName>.pipeline-report.html
```

The report summarizes the configuration, raw read pairs per sample, reads
retained per sample after DADA2, stage-specific read loss, sample-level quality
control, domain composition, 16S/18S/chloroplast/mitochondrial assignments,
unassigned sequences, abundant taxa, and—when enabled—the
internal-standard results and figures.

The taxonomy explorer can filter samples by any populated metadata column in
`config/samples.tsv` and display every taxonomy level present in the formatted
data, from domain through species (including PR2 ranks and ProPortal ecotypes).
It is embedded in the self-contained report and works without an internet
connection.

For convenient downloading, the workflow also copies the report and formatted
data tables into the top-level `Results-Export/` folder. When internal standards
are enabled, the corrected ASV table is included there as well.
