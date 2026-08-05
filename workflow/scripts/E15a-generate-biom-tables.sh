#!/usr/bin/env bash

set -euo pipefail

# Export plaintext taxonomy from both classifiers.
qiime tools export \
    --input-path "${snakemake_input[SILVAclassified]}" \
    --output-path "${snakemake_output[SILVAtaxdir]}"
mv \
    "${snakemake_output[SILVAtaxdir]}/taxonomy.tsv" \
    "${snakemake_output[SILVAtaxfile]}"

qiime tools export \
    --input-path "${snakemake_input[PR2classified]}" \
    --output-path "${snakemake_output[PR2taxdir]}"
mv \
    "${snakemake_output[PR2taxdir]}/taxonomy.tsv" \
    "${snakemake_output[PR2taxfile]}"

sed -i $'1c#OTUID\ttaxonomy\tconfidence' "${snakemake_output[SILVAtaxfile]}"
sed -i $'1c#OTUID\ttaxonomy\tconfidence' "${snakemake_output[PR2taxfile]}"

# Export the complete 18S feature table. Taxonomic subsets are generated from
# the taxonomy-enriched TSVs later so an empty subset remains a valid TSV.
qiime tools export \
    --input-path "${snakemake_input[all18Stable]}" \
    --output-path results/02-euks/15-exports/
mv \
    results/02-euks/15-exports/feature-table.biom \
    "${snakemake_output[all18Stablebiom]}"
