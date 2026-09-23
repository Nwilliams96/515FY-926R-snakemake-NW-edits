#!/usr/bin/env bash

set -euo pipefail

mkdir -p "${snakemake_output[0]}" "$(dirname "${snakemake_log[0]}")"

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${snakemake_input[sequences]}" \
  --p-trim-left 0 \
  --p-trunc-len "${snakemake_params[trunc_len]}" \
  --p-max-ee "${snakemake_params[max_ee]}" \
  --p-trunc-q "${snakemake_params[trunc_q]}" \
  --p-pooling-method "${snakemake_params[pooling_method]}" \
  --p-chimera-method "${snakemake_params[chimera_method]}" \
  --p-min-fold-parent-over-abundance "${snakemake_params[min_fold_parent_over_abundance]}" \
  --p-n-reads-learn "${snakemake_params[n_reads_learn]}" \
  --p-n-threads "${snakemake[threads]}" \
  --o-table "${snakemake_output[table]}" \
  --o-representative-sequences "${snakemake_output[sequences]}" \
  --o-denoising-stats "${snakemake_output[stats]}" \
  --o-base-transition-stats "${snakemake_output[transitions]}" \
  --verbose 2>&1 | tee -a "${snakemake_log[0]}"
