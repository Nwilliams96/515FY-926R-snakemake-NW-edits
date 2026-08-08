#!/usr/bin/env bash

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
  --p-pooling-method ${snakemake_params[pooling_method]} \
  --p-chimera-method ${snakemake_params[chimera_method]} \
  --p-min-fold-parent-over-abundance ${snakemake_params[min_fold_parent_over_abundance]} \
  --p-n-reads-learn ${snakemake_params[n_reads_learn]} \
  --o-table ${snakemake_output[proktable]} \
  --o-representative-sequences ${snakemake_output[prokrepseqs]} \
  --o-denoising-stats ${snakemake_output[prokstats]} \
  --p-n-threads ${snakemake[threads]} \
  --verbose 2>&1 | tee -a ${snakemake_log[0]}
