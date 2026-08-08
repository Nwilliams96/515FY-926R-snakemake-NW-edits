#!/usr/bin/env bash

qiime feature-classifier classify-sklearn \
  --i-classifier ${snakemake_input[classDB]} \
  --i-reads ${snakemake_input[sequences]} \
  --p-n-jobs ${snakemake[threads]} \
  --o-classification ${snakemake_output[classified]}
