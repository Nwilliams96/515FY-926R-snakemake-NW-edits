#!/usr/bin/env bash

set -euo pipefail

biom convert \
    -i "${snakemake_input[all18StablebiomSILVAtax]}" \
    -o "${snakemake_output[all18StablebiomSILVAtaxtsv]}" \
    --to-tsv \
    --header-key taxonomy

biom convert \
    -i "${snakemake_input[all18StablebiomPR2tax]}" \
    -o "${snakemake_output[all18StablebiomPR2taxtsv]}" \
    --to-tsv \
    --header-key taxonomy
