#!/usr/bin/env bash

set -euo pipefail

snakemake \
  --snakefile workflow/Snakefile \
  --use-conda \
  --cores 1 \
  "$@"
