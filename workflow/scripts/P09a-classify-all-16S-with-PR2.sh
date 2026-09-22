#!/usr/bin/env bash
set -euo pipefail

# Classify every unique 16S ASV with PR2.  Confidence is deliberately left at
# zero here so the resolver can record low-confidence plastid candidates while
# applying the user-facing acceptance threshold itself.
qiime feature-classifier classify-sklearn \
  --i-classifier "${snakemake_input[PR2classifier]}" \
  --i-reads "${snakemake_input[prokseqs]}" \
  --p-confidence 0.0 \
  --p-n-jobs "${snakemake[threads]}" \
  --o-classification "${snakemake_output[classified]}" \
  > "${snakemake_log[0]}" 2>&1
