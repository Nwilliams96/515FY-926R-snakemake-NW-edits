#!/usr/bin/env bash

sample=${snakemake_params[0]}
totalPROKseqs=$(gzip -cd "${snakemake_input[prok]}" | awk 'END { print NR / 4 }')
totalEUKseqs=$(gzip -cd "${snakemake_input[euk]}" | awk 'END { print NR / 4 }')

totalSeqs=$(python -c "print($totalEUKseqs + $totalPROKseqs)")
eukFrac=`bc <<< "scale=8; $totalEUKseqs/$totalSeqs"` 

printf "$sample\t$totalPROKseqs\t$totalEUKseqs\t$totalSeqs\t$eukFrac\n" >> ${snakemake_output[eukfrac]}
