#!/usr/bin/env bash

set -euo pipefail

exec bash run_snakemake_USC_CARC_only.sh forward_only_16S_test "$@"
