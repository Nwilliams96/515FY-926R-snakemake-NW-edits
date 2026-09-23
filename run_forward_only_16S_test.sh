#!/usr/bin/env bash

set -euo pipefail

exec bash run_snakemake.sh forward_only_16S_test "$@"
