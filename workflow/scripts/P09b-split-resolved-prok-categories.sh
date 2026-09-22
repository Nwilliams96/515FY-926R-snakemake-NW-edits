#!/usr/bin/env bash
set -euo pipefail

# PR2 plastid lineages contain ":plas"; retained legacy SILVA chloroplast
# lineages contain "o__Chloroplast".  Every split is made from the resolved
# taxonomy so rescued plastids cannot leak into cyanobacterial/prokaryotic bins.
taxonomy="${snakemake_input[resolved_taxonomy]}"
table="${snakemake_input[proktable]}"
plastid_terms=":plas,o__Chloroplast"
cyanobacteria_terms="p__Cyanobacteria,p__Cyanobacteriota"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-include "$plastid_terms" \
  --o-filtered-table "${snakemake_output[includechlorotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-exclude "$plastid_terms" \
  --o-filtered-table "${snakemake_output[excludechlorotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-include "f__Mitochondria" \
  --o-filtered-table "${snakemake_output[onlymitotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-include "$plastid_terms,$cyanobacteria_terms" \
  --o-filtered-table "${snakemake_output[onlyalgaetable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-include "$cyanobacteria_terms" --p-exclude "$plastid_terms" \
  --o-filtered-table "${snakemake_output[onlycyanotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-exclude "f__Mitochondria" \
  --o-filtered-table "${snakemake_output[nomitotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-exclude "f__Mitochondria,$plastid_terms" \
  --o-filtered-table "${snakemake_output[nomitonochlorotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-exclude "f__Mitochondria,$plastid_terms,$cyanobacteria_terms" \
  --o-filtered-table "${snakemake_output[nomitonochloronocyanotable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-include "d__Archaea" \
  --o-filtered-table "${snakemake_output[onlyarchaeatable]}"

qiime taxa filter-table --i-table "$table" --i-taxonomy "$taxonomy" \
  --p-exclude "d__Archaea" \
  --o-filtered-table "${snakemake_output[noarchaeatable]}"
