#!/usr/bin/env bash

set -euo pipefail

trim_overhang_flag="--p-no-trim-overhang"
if [[ "${snakemake_params[trim_overhang]}" == "True" || "${snakemake_params[trim_overhang]}" == "true" ]]; then
  trim_overhang_flag="--p-trim-overhang"
fi

retain_unmerged_flag="--p-no-retain-unmerged"
if [[ "${snakemake_params[retain_unmerged]}" == "True" || "${snakemake_params[retain_unmerged]}" == "true" ]]; then
  retain_unmerged_flag="--p-retain-unmerged"
fi

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${snakemake_input[0]} \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f ${snakemake_params[truncR1]} \
  --p-trunc-len-r ${snakemake_params[truncR2]} \
  --p-max-ee-f ${snakemake_params[max_ee_f]} \
  --p-max-ee-r ${snakemake_params[max_ee_r]} \
  --p-trunc-q ${snakemake_params[trunc_q]} \
  --p-min-overlap ${snakemake_params[min_overlap]} \
  --p-max-merge-mismatch ${snakemake_params[max_merge_mismatch]} \
  ${trim_overhang_flag} \
  ${retain_unmerged_flag} \
  --p-pooling-method ${snakemake_params[pooling_method]} \
  --p-chimera-method ${snakemake_params[chimera_method]} \
  --p-min-fold-parent-over-abundance ${snakemake_params[min_fold_parent_over_abundance]} \
  --p-n-reads-learn ${snakemake_params[n_reads_learn]} \
  --o-table ${snakemake_output[proktable]} \
  --o-representative-sequences ${snakemake_output[prokrepseqs]} \
  --o-denoising-stats ${snakemake_output[prokstats]} \
  --o-base-transition-stats ${snakemake_output[prokbasetransitions]} \
  --p-n-threads ${snakemake[threads]} \
  --verbose 2>&1 | tee -a ${snakemake_log[0]}
