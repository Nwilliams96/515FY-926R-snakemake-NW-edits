#!/usr/bin/env bash

classifier="${snakemake_input[classDB]}"
if [[ -n "${snakemake_params[classDB]}" ]]; then
  classifier="${snakemake_params[classDB]}"
fi

qiime feature-classifier classify-sklearn \
  --i-classifier "${classifier}" \
  --i-reads ${snakemake_input[sequences]} \
  --p-n-jobs ${snakemake[threads]} \
  --o-classification ${snakemake_output[classified]}
