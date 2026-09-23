"""Estimate primer-to-primer lengths for forward-only-unique ASVs.

The estimates are reference based.  Each ASV is aligned to SILVA's full-length
DADA2 species reference, then the configured primer sites are located around
the aligned query coordinates.  Results must not be interpreted as direct
measurements of the environmental molecules.
"""

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


IUPAC = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "U": frozenset("T"), "R": frozenset("AG"),
    "Y": frozenset("CT"), "S": frozenset("GC"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"), "B": frozenset("CGT"),
    "D": frozenset("AGT"), "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVN", "TGCAYRKMVHDBN"
)
BLAST_COLUMNS = [
    "query", "subject", "identity", "alignment_length", "query_length",
    "mismatches", "gap_opens", "query_start", "query_end",
    "subject_start", "subject_end", "evalue", "bitscore",
]


def reverse_complement(sequence):
    return sequence.upper().translate(COMPLEMENT)[::-1]


def mismatch_count(reference, motif):
    mismatches = 0
    for observed, expected in zip(reference.upper(), motif.upper()):
        if not IUPAC.get(observed, frozenset(observed)) & IUPAC[expected]:
            mismatches += 1
    return mismatches


def best_site(sequence, motif, start, end, expected_position):
    motif_length = len(motif)
    start = max(0, int(start))
    end = min(len(sequence) - motif_length + 1, int(end))
    if end <= start:
        return None
    candidates = []
    for position in range(start, end):
        mismatches = mismatch_count(
            sequence[position:position + motif_length], motif
        )
        candidates.append((mismatches, abs(position - expected_position), position))
    return min(candidates) if candidates else None


def read_blast(path):
    hits = defaultdict(list)
    with Path(path).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=BLAST_COLUMNS)
        for row in reader:
            if not row.get("query"):
                continue
            for key in (
                "identity", "evalue", "bitscore",
            ):
                row[key] = float(row[key])
            for key in (
                "alignment_length", "query_length", "mismatches", "gap_opens",
                "query_start", "query_end", "subject_start", "subject_end",
            ):
                row[key] = int(row[key])
            hits[row["query"]].append(row)
    for query in hits:
        hits[query].sort(key=lambda row: (-row["bitscore"], -row["identity"], row["evalue"]))
    return hits


def read_abundances(path):
    abundances = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Forward_only_unique_ASV", "").lower() == "true":
                abundances[row["Forward_ASV_ID"]] = float(row["Total_forward_reads"])
    return abundances


def load_references(path, wanted):
    records = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        identifier = None
        header = ""
        sequence = []
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if identifier in wanted:
                    records[identifier] = (header, "".join(sequence).upper().replace("U", "T"))
                header = line[1:]
                identifier = header.split()[0]
                sequence = []
            elif identifier in wanted:
                sequence.append(line)
        if identifier in wanted:
            records[identifier] = (header, "".join(sequence).upper().replace("U", "T"))
    return records


def oriented_reference(sequence, hit):
    if hit["subject_start"] <= hit["subject_end"]:
        # BLAST coordinates are one-based; this is the reference coordinate at
        # which query position query_start begins.
        aligned_start = hit["subject_start"] - 1
        return sequence, aligned_start
    aligned_start = len(sequence) - hit["subject_start"]
    return reverse_complement(sequence), aligned_start


def evaluate_hit(hit, reference, forward_primer, reverse_binding_motif,
                 trunc_f, trunc_r, min_overlap):
    header, raw_sequence = reference
    sequence, aligned_start = oriented_reference(raw_sequence, hit)
    query_reference_start = aligned_start - (hit["query_start"] - 1)

    expected_forward = query_reference_start - len(forward_primer)
    forward = best_site(
        sequence, forward_primer,
        expected_forward - 40, expected_forward + 41,
        expected_forward,
    )
    expected_reverse = query_reference_start + 370
    reverse = best_site(
        sequence, reverse_binding_motif,
        query_reference_start + 250, query_reference_start + 651,
        expected_reverse,
    )
    if forward is None or reverse is None:
        return None

    f_mismatch, _, f_start = forward
    r_mismatch, _, r_start = reverse
    max_f_mismatch = math.floor(len(forward_primer) * 0.2)
    max_r_mismatch = math.floor(len(reverse_binding_motif) * 0.2)
    if f_mismatch > max_f_mismatch or r_mismatch > max_r_mismatch:
        return None

    post_primer_length = r_start - (f_start + len(forward_primer))
    if not 100 <= post_primer_length <= 1000:
        return None
    amplicon_length = (
        r_start + len(reverse_binding_motif) - f_start
    )
    predicted_overlap = trunc_f + trunc_r - post_primer_length
    return {
        "Reference_ID": hit["subject"],
        "Reference_header": header,
        "BLAST_identity_percent": f"{hit['identity']:.3f}",
        "BLAST_query_coverage_percent": f"{100 * hit['alignment_length'] / hit['query_length']:.3f}",
        "BLAST_evalue": f"{hit['evalue']:.6g}",
        "Forward_primer_mismatches": f_mismatch,
        "Reverse_primer_mismatches": r_mismatch,
        "Estimated_amplicon_length_including_primers": amplicon_length,
        "Estimated_post_primer_insert_length": post_primer_length,
        "Predicted_overlap_at_configured_truncation": predicted_overlap,
        "Length_mergeable_at_configured_min_overlap": str(predicted_overlap >= min_overlap).lower(),
        "Reference_estimate_status": "estimated",
        "_score": (f_mismatch + r_mismatch, -hit["bitscore"], -hit["identity"]),
    }


def write_tsv(path, columns, rows):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def run(snakemake_object):
    hits = read_blast(snakemake_object.input.hits)
    abundances = read_abundances(snakemake_object.input.comparison)
    wanted_references = {
        hit["subject"] for query_hits in hits.values() for hit in query_hits
    }
    references = load_references(snakemake_object.input.reference, wanted_references)

    forward_primer = str(snakemake_object.params.forward_primer).upper().replace("U", "T")
    reverse_primer = str(snakemake_object.params.reverse_primer).upper().replace("U", "T")
    reverse_binding_motif = reverse_complement(reverse_primer)
    trunc_f = int(snakemake_object.params.trunc_len_f)
    trunc_r = int(snakemake_object.params.trunc_len_r)
    min_overlap = int(snakemake_object.params.min_overlap)

    rows = []
    for query in sorted(abundances, key=lambda value: (-abundances[value], value)):
        candidates = []
        for hit in hits.get(query, []):
            reference = references.get(hit["subject"])
            if reference is None:
                continue
            estimate = evaluate_hit(
                hit, reference, forward_primer, reverse_binding_motif,
                trunc_f, trunc_r, min_overlap,
            )
            if estimate is not None:
                candidates.append(estimate)

        if candidates:
            selected = min(candidates, key=lambda row: row["_score"])
        else:
            selected = {
                "Reference_estimate_status": (
                    "no_qualifying_reference_hit" if not hits.get(query)
                    else "primer_sites_not_resolved"
                )
            }
        selected = dict(selected)
        selected.pop("_score", None)
        rows.append(
            {
                "Forward_ASV_ID": query,
                "Total_forward_reads": (
                    str(int(abundances[query])) if abundances[query].is_integer()
                    else f"{abundances[query]:.10g}"
                ),
                **selected,
            }
        )

    columns = [
        "Forward_ASV_ID", "Total_forward_reads", "Reference_ID",
        "Reference_header", "BLAST_identity_percent",
        "BLAST_query_coverage_percent", "BLAST_evalue",
        "Forward_primer_mismatches", "Reverse_primer_mismatches",
        "Estimated_amplicon_length_including_primers",
        "Estimated_post_primer_insert_length",
        "Predicted_overlap_at_configured_truncation",
        "Length_mergeable_at_configured_min_overlap",
        "Reference_estimate_status",
    ]
    write_tsv(snakemake_object.output.estimates, columns, rows)

    estimated = [row for row in rows if row.get("Reference_estimate_status") == "estimated"]
    too_long = [
        row for row in estimated
        if row.get("Length_mergeable_at_configured_min_overlap") == "false"
    ]
    total_unique_reads = sum(abundances.values())
    too_long_reads = sum(float(row["Total_forward_reads"]) for row in too_long)
    summary = [
        {"Metric": "forward_only_unique_ASVs", "Value": len(rows)},
        {"Metric": "forward_only_unique_reads", "Value": f"{total_unique_reads:.10g}"},
        {"Metric": "ASVs_with_reference_length_estimate", "Value": len(estimated)},
        {"Metric": "ASVs_predicted_too_long_to_merge", "Value": len(too_long)},
        {"Metric": "reads_in_ASVs_predicted_too_long_to_merge", "Value": f"{too_long_reads:.10g}"},
        {"Metric": "maximum_mergeable_post_primer_length", "Value": trunc_f + trunc_r - min_overlap},
        {"Metric": "interpretation", "Value": "Reference-based estimates are hypotheses, not direct insert-length measurements."},
    ]
    write_tsv(snakemake_object.output.summary, ["Metric", "Value"], summary)


if "snakemake" in globals():
    run(snakemake)
