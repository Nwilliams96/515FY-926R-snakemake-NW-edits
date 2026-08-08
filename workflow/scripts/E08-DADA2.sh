#!/usr/bin/env bash

#230904 added in --p-trunc-q 0
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ${snakemake_input[0]} \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --p-max-ee ${snakemake_params[max_ee]} \
  --p-trunc-q ${snakemake_params[trunc_q]} \
  --p-pooling-method ${snakemake_params[pooling_method]} \
  --p-chimera-method ${snakemake_params[chimera_method]} \
  --p-min-fold-parent-over-abundance ${snakemake_params[min_fold_parent_over_abundance]} \
  --p-n-reads-learn ${snakemake_params[n_reads_learn]} \
  --p-n-threads ${snakemake[threads]} \
  --o-table ${snakemake_output[euktable]} \
  --o-representative-sequences ${snakemake_output[eukrepseqs]} \
  --o-denoising-stats ${snakemake_output[eukstats]} \
  --verbose
