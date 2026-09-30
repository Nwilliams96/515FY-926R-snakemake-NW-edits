"""Create one FASTA file per configured internal standard."""

import csv
import re
import warnings
from pathlib import Path


IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


def reverse_complement(sequence):
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


def find_iupac(sequence, pattern, start=0):
    """Return the first ambiguity-aware pattern match at or after start."""
    last_start = len(sequence) - len(pattern)
    for offset in range(start, last_start + 1):
        if all(
            IUPAC_BASES[sequence[offset + index]] & IUPAC_BASES[base]
            for index, base in enumerate(pattern)
        ):
            return offset
    return None


def primer_bounded_insert(sequence, forward_primer, reverse_primer):
    """Return the insert in forward orientation, excluding both primers."""
    reverse_site = reverse_complement(reverse_primer)
    for oriented_sequence in (sequence, reverse_complement(sequence)):
        forward_start = find_iupac(oriented_sequence, forward_primer)
        if forward_start is None:
            continue
        insert_start = forward_start + len(forward_primer)
        reverse_start = find_iupac(oriented_sequence, reverse_site, insert_start)
        if reverse_start is not None:
            return oriented_sequence[insert_start:reverse_start]
    return None


forward_primer = str(snakemake.params.forward_primer).strip().upper()
reverse_primer = str(snakemake.params.reverse_primer).strip().upper()
trunc_r1 = int(snakemake.params.trunc_r1)
trunc_r2 = int(snakemake.params.trunc_r2)
if trunc_r1 <= 0 or trunc_r2 <= 0:
    raise ValueError("18S concatenation lengths must both be positive")


with open(snakemake.input[0], newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    sequence_column = next(
        (
            column
            for column in ("full_SSU_sequence", "full_16S_sequence")
            if column in (reader.fieldnames or [])
        ),
        None,
    )
    required_columns = {"internal_std_ID"}
    missing_columns = required_columns.difference(reader.fieldnames or [])
    if sequence_column is None:
        missing_columns.add("full_SSU_sequence (or legacy full_16S_sequence)")
    if missing_columns:
        raise ValueError(
            "config/internal_stds.tsv is missing column(s): "
            + ", ".join(sorted(missing_columns))
        )
    standard_rows = list(reader)
    standards = {row["internal_std_ID"].strip(): row for row in standard_rows}
    if len(standards) != len(standard_rows):
        raise ValueError("config/internal_stds.tsv contains duplicate internal_std_ID values")


configured_ids = [
    str(standard_id).strip() for standard_id in snakemake.params.standard_ids
]
if not configured_ids:
    raise ValueError("At least one internal standard must be configured")
if len(set(configured_ids)) != len(configured_ids):
    raise ValueError("Configured internal standard names must be unique")
for standard_id in configured_ids:
    if not re.fullmatch(r"[A-Za-z0-9._-]+", standard_id):
        raise ValueError(
            f"Internal standard name {standard_id!r} may contain only letters, "
            "numbers, periods, underscores, and hyphens"
        )


output_paths = list(snakemake.output.fastas)
if len(output_paths) != len(configured_ids):
    raise ValueError("Expected one FASTA output per configured internal standard")

for standard_id, output_name in zip(configured_ids, output_paths):
    if standard_id not in standards:
        raise ValueError(
            f"Configured internal standard {standard_id!r} is not present in "
            "config/internal_stds.tsv"
        )

    sequence = re.sub(
        r"\s+", "", standards[standard_id][sequence_column]
    ).upper()
    if not sequence:
        raise ValueError(f"Internal standard {standard_id!r} has an empty sequence")
    if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
        raise ValueError(
            f"Internal standard {standard_id!r} contains invalid nucleotide symbols"
        )

    output_path = Path(output_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    records = [f">{standard_id}\n{sequence}\n"]

    # The 18S workflow does not merge overlapping reads. It trims R1 and R2,
    # reverse-complements R2, and joins the two ends directly with no spacer.
    # Build the same artificial sequence so a denoised 18S ISD ASV can obtain
    # 100% query coverage in BLAST despite the omitted middle of the amplicon.
    insert = primer_bounded_insert(sequence, forward_primer, reverse_primer)
    if insert is None:
        warnings.warn(
            f"Could not find an ordered {forward_primer}/{reverse_primer} primer "
            f"pair in internal standard {standard_id!r}; only its continuous "
            "reference will be available for BLAST"
        )
    elif len(insert) < max(trunc_r1, trunc_r2):
        warnings.warn(
            f"The primer-bounded insert for internal standard {standard_id!r} "
            "is shorter than an 18S truncation length; only its continuous "
            "reference will be available for BLAST"
        )
    else:
        concatenated = insert[:trunc_r1] + insert[-trunc_r2:]
        records.append(
            f">{standard_id}__18S_concatenated_{trunc_r1}_{trunc_r2}\n"
            f"{concatenated}\n"
        )

    output_path.write_text("".join(records), encoding="utf-8")
