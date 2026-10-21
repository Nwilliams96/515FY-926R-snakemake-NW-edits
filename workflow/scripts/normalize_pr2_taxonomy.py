"""Remove PR2 rank prefixes without changing confidence or other TSV columns.

Downloaded and existing classifiers can retain RESCRIPt-style rank prefixes.
The established downstream tables use nine unprefixed PR2 ranks and identify
SILVA lineages by their d__ prefix. Normalize only taxonomy text after export.
RESCRIPt edit-taxonomy is suitable for reference lineages but returns a Series,
so using it on classification results would discard their Confidence column.

"""

import csv
import re
from pathlib import Path


# Permit spacing used by externally trained classifiers while preserving it.
# A prefix-like substring inside a name is not a rank marker and stays intact.
RANK_PREFIX = re.compile(r"(^|;)(\s*)(?:d__|sgr__|dv__|dvs__|c__|o__|f__|g__|s__)")


def normalize_pr2_taxonomy(input_path, output_path):
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    changed = 0
    with Path(input_path).open(encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError("The exported PR2 taxonomy TSV is empty")
        taxonomy_columns = [
            index for index, name in enumerate(header)
            if name.strip().lower() in ("taxon", "taxonomy")
        ]
        if len(taxonomy_columns) != 1:
            raise ValueError("Expected one Taxon or taxonomy column in the PR2 TSV")
        column = taxonomy_columns[0]
        with output.open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for row_number, row in enumerate(reader, 2):
                if not row:
                    writer.writerow(row)
                    continue
                if len(row) != len(header):
                    raise ValueError(f"PR2 taxonomy row {row_number} has an unexpected number of columns")
                normalized = RANK_PREFIX.sub(r"\1\2", row[column])
                changed += normalized != row[column]
                row[column] = normalized
                # csv operates on strings: IDs such as 001 and Confidence text
                # such as 0.8000 are preserved, not parsed and rounded.
                writer.writerow(row)
    return changed


if "snakemake" in globals():
    changed = normalize_pr2_taxonomy(
        snakemake.input.taxonomy, snakemake.output.normalized
    )
    if snakemake.log:
        log = Path(snakemake.log[0])
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            f"Removed PR2 rank prefixes from {changed} taxonomy row(s); all other columns retained.\n",
            encoding="utf-8",
        )
