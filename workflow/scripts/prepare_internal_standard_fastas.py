"""Create one FASTA file per configured internal standard."""

import csv
import re
from pathlib import Path


with open(snakemake.input[0], newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required_columns = {"internal_std_ID", "full_16S_sequence"}
    missing_columns = required_columns.difference(reader.fieldnames or [])
    if missing_columns:
        raise ValueError(
            "config/internal_stds.tsv is missing column(s): "
            + ", ".join(sorted(missing_columns))
        )
    standards = {row["internal_std_ID"].strip(): row for row in reader}


for slot in ("intstd1", "intstd2", "intstd3"):
    standard_id = str(snakemake.params[f"{slot}name"]).strip()
    if standard_id not in standards:
        raise ValueError(
            f"Internal standard {standard_id!r} is configured as {slot} but is "
            "not present in config/internal_stds.tsv"
        )

    sequence = re.sub(
        r"\s+", "", standards[standard_id]["full_16S_sequence"]
    ).upper()
    if not sequence:
        raise ValueError(f"Internal standard {standard_id!r} has an empty sequence")
    if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
        raise ValueError(
            f"Internal standard {standard_id!r} contains invalid nucleotide symbols"
        )

    output_path = Path(snakemake.output[slot])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(f">{standard_id}\n{sequence}\n", encoding="utf-8")
