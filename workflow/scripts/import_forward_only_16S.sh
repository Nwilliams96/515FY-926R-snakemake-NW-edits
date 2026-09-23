#!/usr/bin/env bash

set -euo pipefail

mkdir -p "$(dirname "${snakemake_output[sequences]}")" "$(dirname "${snakemake_log[0]}")"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${snakemake_input[manifest]}" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path "${snakemake_output[sequences]}" \
  2>&1 | tee -a "${snakemake_log[0]}"
