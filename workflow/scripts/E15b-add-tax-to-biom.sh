#!/usr/bin/env bash

set -euo pipefail

biom add-metadata \
    -i "${snakemake_input[all18Stablebiom]}" \
    -o "${snakemake_output[all18StablebiomPR2tax]}" \
    --observation-metadata-fp "${snakemake_input[PR2taxfile]}" \
    --sc-separated taxonomy
