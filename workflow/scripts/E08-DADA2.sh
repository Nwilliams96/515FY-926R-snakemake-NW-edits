#!/usr/bin/env bash

# AI-assisted change — see provenance ledger linked from README.md, commits f7aba4f and 5c6c1aa.
# Purpose: pass configured single-end 18S DADA2 controls and retain transition
# statistics needed by QIIME 2 2026.7 and the final report.
# Human review status: pending.

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
  --o-base-transition-stats ${snakemake_output[eukbasetransitions]} \
  --verbose
