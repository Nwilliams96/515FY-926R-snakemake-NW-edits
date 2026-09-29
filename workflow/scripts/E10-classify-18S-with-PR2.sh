#!/usr/bin/env bash

set -euo pipefail

mkdir -p "$(dirname "${snakemake_output[classified]}")"

# Classify 18S ASVs once with PR2.
qiime feature-classifier classify-sklearn \
  --i-classifier "${snakemake_params[classDB]}" \
  --i-reads "${snakemake_input[sequences]}" \
  --p-n-jobs "${snakemake[threads]}" \
  --o-classification "${snakemake_output[classified]}"

# Preserve the established workaround for whitespace in PR2 taxonomy metadata.
qiime tools export \
  --input-path "${snakemake_output[classified]}" \
  --output-path "${snakemake_output[exported_taxonomy]}"

qiime metadata tabulate \
  --m-input-file "${snakemake_output[exported_taxonomy]}/taxonomy.tsv" \
  --o-visualization "${snakemake_output[taxonomy_metadata]}"

qiime tools export \
  --input-path "${snakemake_output[taxonomy_metadata]}" \
  --output-path "${snakemake_output[exported_metadata]}"

sed -i '/#q2:types/d' "${snakemake_output[exported_metadata]}/metadata.tsv"

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path "${snakemake_output[exported_metadata]}/metadata.tsv" \
  --output-path "${snakemake_output[normalized]}"
