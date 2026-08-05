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
    standard_rows = list(reader)
    standards = {row["internal_std_ID"].strip(): row for row in standard_rows}
    if len(standards) != len(standard_rows):
        raise ValueError("config/internal_stds.tsv contains duplicate internal_std_ID values")


configured_ids = {
    slot: str(snakemake.params[f"{slot}name"]).strip()
    for slot in ("intstd1", "intstd2", "intstd3")
}
if len(set(configured_ids.values())) != 3:
    raise ValueError("The three configured internal standard names must be unique")
for standard_id in configured_ids.values():
    if not re.fullmatch(r"[A-Za-z0-9._-]+", standard_id):
        raise ValueError(
            f"Internal standard name {standard_id!r} may contain only letters, "
            "numbers, periods, underscores, and hyphens"
        )


for slot, standard_id in configured_ids.items():
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
