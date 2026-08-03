#!/usr/bin/env bash

makeblastdb -dbtype nucl -in ${snakemake_input[intstd1]}
touch ${snakemake_output[intstd1]}

makeblastdb -dbtype nucl -in ${snakemake_input[intstd2]}
touch ${snakemake_output[intstd2]}

makeblastdb -dbtype nucl -in ${snakemake_input[intstd3]}
touch ${snakemake_output[intstd3]}
