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

The required starting-pool molarity file is named
`config/prok_and_euk_SSU_amplicon_molarities.tsv`. Updated workflows still
accept the previous `config/prok_and_euk_SSU_amplicon_concentrations.tsv` or
`config/bioanalyzer.tsv` when the new file is absent, allowing older study
folders to finish after a pipeline update.

Additional user-defined columns may be added to `config/samples.tsv`. They are
preserved in the generated QIIME metadata and are available as filters in the
HTML report; they do not alter read processing.

The cloned repository intentionally does not include a `config/` folder. Create
and download the study-specific setup package from the pipeline tutorial, then
place its complete `config/` folder in the repository before running Snakemake.

## Taxonomy databases

Prokaryotic ASVs are classified with the official SILVA 144 SSU Ref NR99
uniform QIIME 2 classifier. With the pipeline's default 515Y/926R primers, the
workflow downloads SILVA's matching V4–V5 region classifier and verifies its
published checksum before placing it in the shared `database_dir`. The bundled
QIIME 2 2026.7 environment matches the version used to train that classifier.
Older generated configs that still name this repository's bundled QIIME 2
2024.5 or 2025.7 definition are automatically migrated to the 2026.7 definition.
The bundled Linux definition is based on the official QIIME 2 2026.7 file but
omits Deblur and its obsolete SortMeRNA 2.0 dependency. This workflow uses
DADA2 exclusively, and excluding that unused dependency prevents Conda solver
failures without changing any pipeline analysis step.

SILVA 144 introduces a prokaryotic `Kingdom` rank and provides a consistent
seven-rank lineage from Domain through Genus. The formatted long table and HTML
report therefore include `Kingdom`. Taxonomy parsing is based on rank prefixes
(`d__`, `k__`, `p__`, and so on), so that new rank cannot shift Phylum, Class,
Order, Family, or Genus into the wrong columns. The official uniform QIIME 2
classifier stops at Genus. The workflow therefore follows classification with
an exact-match search against SILVA's official version 144 DADA2 species
reference. An unambiguous exact match is accepted only when its genus agrees
with the QIIME 2 assignment. For example, an accepted match is appended as
`s__Vibrio cholerae`, and the formatted outputs contain `Vibrio` in `Genus` and
`Vibrio cholerae` in `Species`. Unmatched, ambiguous, or genus-conflicting ASVs
remain blank at Species. SILVA notes that these organism names are not curated,
so species calls should be treated as provisional. PR2-derived eukaryotic
species labels remain available independently.

The chloroplast/cyanobacterial subsetting step accepts both the older
`p__Cyanobacteria` label and SILVA 144's `p__Cyanobacteriota` label, preserving
those downstream tables across the database transition.

If non-default primers are entered, the workflow uses SILVA's official
full-length classifier instead of applying the 515Y/926R-specific model to an
incompatible region. PR2 remains the primary eukaryotic classifier.

Both release-pinned classifiers and the SILVA species reference live in the
shared `database_dir`. After this upgrade, the first run downloads the SILVA
144 classifier, downloads the approximately 141 MB species reference, and
builds a QIIME 2 2026.7-labelled PR2 classifier if they are not already present;
later project clones reuse them. This one-time database setup also occurs when
the corresponding pre-existing setting is `false`; later project clones reuse
them. Database families are controlled independently with
`use_preexisting_bbsplit_database`, `use_preexisting_silva_database`, and
`use_preexisting_pr2_database`. Marking one as `true` requires its release- and
primer-compatible files to exist under `database_dir`; marking it as `false`
allows Snakemake to prepare it if missing. The older
`use_preexisting_databases` setting remains supported as a BBsplit fallback.

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

The small per-sample BBTools trimming, repair, and fusion commands use a fixed
2 GB Java heap and declare 2.5 GB per job to Snakemake. This prevents concurrent
BBTools processes from each auto-claiming most of the node's available memory.

Sample identifiers are imported as text when the final QIIME metadata tables
are assembled, so numeric-only IDs and leading zeroes are retained correctly.
DADA2 export scripts also propagate failures from QIIME and BIOM commands rather
than allowing a failed piped command to appear successful.

When databases must be built, the official pre-trained SILVA classifier is
downloaded while PR2 preparation proceeds independently. BBSplit indexing and
PR2 primer extraction use up to eight cores. QIIME's naive-Bayes classifier-
training action does not expose a worker-count option, so the final PR2
classifier-training command remains single-core.

## Run report

The final workflow target creates a self-contained report at:

```text
results/07-report/<studyName>.pipeline-report.html
```

The report summarizes the configuration, raw read pairs per sample, primer
trimming and BBsplit assignment loss, reads retained per sample after DADA2,
stage-specific DADA2 loss, sample-level quality control, domain composition,
16S/18S/chloroplast/mitochondrial assignments, unassigned sequences, abundant
taxa, the exact 16S and 18S molarity/read correction factors used during
merging, and—when enabled—the internal-standard results and figures. Horizontally
scrollable sample figures keep their y-axis scale visible, and report palettes
are chosen for colour-vision accessibility.

The taxonomy explorer provides a QIIME 2-style, 100%-stacked taxonomic bar plot
with one bar per sample. It can order and filter plotted samples by SampleID,
Condition, Latitude, Longitude, or Depth from `config/samples.tsv`, and display
every taxonomy level present in the formatted data, from Domain through Species,
including SILVA 144's Kingdom rank, PR2 ranks, and ProPortal ecotypes. It is embedded in the
self-contained report and works without an internet connection.

For convenient downloading, the workflow copies the report, formatted data,
and phylum- and order-level summaries into a top-level
`<projectName>-Results-Export/` folder. When internal standards are enabled, the
corrected ASV table is included there as well. Older configs without
`projectName` use `studyName` for the folder prefix.

When internal standards are enabled, the ordinary exported `long_data.tsv` is
also rebuilt without every ASV identified as an internal standard. The report
uses this same filtered table, preventing ISD reads from appearing in domain or
taxonomic-abundance plots. The separate `ISD_corrected_long_data.tsv` retains
the absolute-copy correction columns but likewise excludes the standard ASVs.
