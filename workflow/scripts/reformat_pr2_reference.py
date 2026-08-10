"""Give PR2 sequences compact feature IDs and export their taxonomy headers."""

from pathlib import Path


def reformat_pr2_reference(input_path, clean_path, taxonomy_path):
    feature_number = 0
    with open(input_path, encoding="utf-8") as source, open(
        clean_path, "w", encoding="utf-8"
    ) as clean, open(taxonomy_path, "w", encoding="utf-8") as taxonomy:
        for line in source:
            if line.startswith(">"):
                feature_number += 1
                original_header = line[1:].strip()
                feature_id = f"feature_{feature_number}"
                clean.write(f">{feature_id}\n")
                taxonomy.write(f"{feature_id}\t{original_header}\n")
            else:
                clean.write(line)
    if feature_number == 0:
        raise ValueError("The PR2 reference FASTA contains no sequence headers")


def run_from_snakemake(snakemake_object):
    reformat_pr2_reference(
        str(snakemake_object.input[0]),
        str(snakemake_object.output.clean),
        str(snakemake_object.output.headers),
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
