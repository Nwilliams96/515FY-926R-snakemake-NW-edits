#!/usr/bin/env bash

qiime feature-classifier classify-sklearn \
  --i-classifier ${snakemake_input[classDB]} \
  --i-reads ${snakemake_input[sequences]} \
  --o-classification ${snakemake_output[classified]}
