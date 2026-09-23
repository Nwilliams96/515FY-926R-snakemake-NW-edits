#!/usr/bin/env bash

set -euo pipefail

mkdir -p "$(dirname "${snakemake_output[hits]}")" "$(dirname "${snakemake_log[0]}")"

if ! grep -q '^>' "${snakemake_input[query]}"; then
  : > "${snakemake_output[hits]}"
  printf 'No forward-only-unique ASVs were available for reference alignment.\n' > "${snakemake_log[0]}"
  exit 0
fi

work_dir="$(mktemp -d)"
trap 'rm -rf "${work_dir}"' EXIT
reference_fasta="${work_dir}/silva-species-reference.fasta"

gzip -dc "${snakemake_input[reference]}" > "${reference_fasta}"
makeblastdb -dbtype nucl -in "${reference_fasta}" \
  > "${snakemake_log[0]}" 2>&1

blastn \
  -task megablast \
  -query "${snakemake_input[query]}" \
  -db "${reference_fasta}" \
  -perc_identity 90 \
  -qcov_hsp_perc 80 \
  -max_target_seqs 10 \
  -num_threads "${snakemake[threads]}" \
  -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore' \
  -out "${snakemake_output[hits]}" \
  2>> "${snakemake_log[0]}"
