"""Create a QIIME 2 single-end manifest from existing prokaryotic R1 FASTQs."""

import csv
from pathlib import Path


def run(snakemake_object):
    output = Path(str(snakemake_object.output[0]))
    output.parent.mkdir(parents=True, exist_ok=True)
    paths = [Path(str(path)).resolve() for path in snakemake_object.input.r1]
    sample_ids = [str(value) for value in snakemake_object.params.sample_ids]

    if len(paths) != len(sample_ids):
        raise ValueError("Forward FASTQ and sample-ID counts do not match")

    rows = []
    for sample_id, path in zip(sample_ids, paths):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty forward FASTQ: {path}")
        rows.append((sample_id.replace("_", "-"), str(path)))

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample-id", "absolute-filepath"])
        writer.writerows(rows)


if "snakemake" in globals():
    run(snakemake)
