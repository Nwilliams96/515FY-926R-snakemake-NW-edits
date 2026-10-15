"""Export linked 16S features retained only by the experimental DADA2 run."""

import csv
import re
from collections import Counter
from pathlib import Path


# AI-assisted change — see CODE_PROVENANCE.md, commit f13649a.
# Purpose: construct the optional linked/unmerged 16S supplement after the
# QIIME and BIOM exports have been made visible as individual Snakemake rules.
# Human review status: pending.


AUDIT_COLUMNS = [
    "ASV_hash", "SILVA_taxonomy", "SILVA_confidence", "SILVA_chloroplast",
    "PR2_taxonomy", "PR2_confidence", "PR2_plastid",
    "PR2_override_accepted", "Final_taxonomy", "Final_confidence",
    "Chloroplast_detection_source",
]


def first_value(row, names):
    for name in names:
        if name in row:
            return (row.get(name) or "").strip()
    return ""


def confidence(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_taxonomy(path):
    rows = []
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            feature_id = first_value(row, ["Feature ID", "FeatureID", "feature-id", "ID"])
            if feature_id and not feature_id.startswith("#"):
                confidence_text = first_value(row, ["Confidence", "confidence"])
                rows.append(
                    {
                        "feature_id": feature_id,
                        "taxonomy": first_value(row, ["Taxon", "Taxonomy", "taxon"]),
                        "confidence_text": confidence_text,
                        "confidence": confidence(confidence_text),
                    }
                )
    return rows


def resolve_taxonomies(silva_rows, pr2_rows, min_confidence):
    pr2_by_id = {row["feature_id"]: row for row in pr2_rows}
    audit_rows = []
    counts = Counter()
    for silva in silva_rows:
        pr2 = pr2_by_id.get(
            silva["feature_id"],
            {"taxonomy": "", "confidence_text": "", "confidence": None},
        )
        silva_chloroplast = bool(
            re.search(r"(?:^|;\s*)o__Chloroplast(?:\s*;|\s*$)", silva["taxonomy"], re.I)
        )
        pr2_plastid = ":plas" in pr2["taxonomy"].lower()
        accepted = bool(
            pr2_plastid
            and pr2["confidence"] is not None
            and pr2["confidence"] >= min_confidence
        )
        if accepted:
            final_taxonomy = pr2["taxonomy"]
            final_confidence = pr2["confidence_text"]
            source = "PR2_confirmed_SILVA" if silva_chloroplast else "PR2_rescue"
        else:
            final_taxonomy = silva["taxonomy"]
            final_confidence = silva["confidence_text"]
            if silva_chloroplast:
                source = "SILVA_chloroplast_unconfirmed"
            elif pr2_plastid:
                source = "PR2_low_confidence_not_accepted"
            else:
                source = "SILVA_retained"
        audit_rows.append(
            {
                "ASV_hash": silva["feature_id"],
                "SILVA_taxonomy": silva["taxonomy"],
                "SILVA_confidence": silva["confidence_text"],
                "SILVA_chloroplast": str(silva_chloroplast).lower(),
                "PR2_taxonomy": pr2["taxonomy"],
                "PR2_confidence": pr2["confidence_text"],
                "PR2_plastid": str(pr2_plastid).lower(),
                "PR2_override_accepted": str(accepted).lower(),
                "Final_taxonomy": final_taxonomy,
                "Final_confidence": final_confidence,
                "Chloroplast_detection_source": source,
            }
        )
        counts[source] += 1
    return audit_rows


def read_biom_tsv(path):
    """Return feature-by-sample counts from a BIOM-converted TSV."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if not line.startswith("# Constructed")]
    reader = csv.DictReader(lines, delimiter="\t")
    rows = {}
    for row in reader:
        feature_id = (row.pop("#OTU ID", "") or row.pop("Feature ID", "")).strip()
        if not feature_id:
            continue
        rows[feature_id] = {
            sample: float(value or 0)
            for sample, value in row.items()
            if sample != "taxonomy"
        }
    return rows


def read_fasta(path):
    records = {}
    identifier = None
    pieces = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(pieces)
                identifier = line[1:].split()[0]
                pieces = []
            elif identifier is not None:
                pieces.append(line)
    if identifier is not None:
        records[identifier] = "".join(pieces)
    return records


def build_supplement(retained_counts, standard_counts, sequences, silva_rows, pr2_rows, min_confidence):
    """Build a wide table and one-row recovery summary."""
    rescued_ids = sorted(set(retained_counts) - set(standard_counts))
    silva_rescued = [row for row in silva_rows if row["feature_id"] in rescued_ids]
    pr2_rescued = [row for row in pr2_rows if row["feature_id"] in rescued_ids]
    audit_rows = resolve_taxonomies(silva_rescued, pr2_rescued, min_confidence)
    audit_by_id = {row["ASV_hash"]: row for row in audit_rows}
    samples = sorted(
        {sample for feature_id in rescued_ids for sample in retained_counts[feature_id]}
    )
    rows = []
    for feature_id in rescued_ids:
        audit = audit_by_id.get(feature_id, {})
        row = {
            "Feature_ID": feature_id,
            "Sequence_Representation": "unmerged_linked_pair",
            "Linked_Sequence": sequences.get(feature_id, ""),
        }
        row.update({column: audit.get(column, "") for column in AUDIT_COLUMNS[1:]})
        row.update(
            {
                sample: f"{retained_counts[feature_id].get(sample, 0):g}"
                for sample in samples
            }
        )
        rows.append(row)

    rescued_reads = sum(
        sum(retained_counts[feature_id].values()) for feature_id in rescued_ids
    )
    chloroplast_ids = {
        row["ASV_hash"]
        for row in audit_rows
        if row.get("Chloroplast_detection_source")
        in {"PR2_confirmed_SILVA", "PR2_rescue", "SILVA_chloroplast_unconfirmed"}
    }
    summary = {
        "standard_non_chimeric_16S_reads": f"{sum(sum(row.values()) for row in standard_counts.values()):g}",
        "retained_run_non_chimeric_16S_reads": f"{sum(sum(row.values()) for row in retained_counts.values()):g}",
        "supplemental_unmerged_16S_features": len(rescued_ids),
        "supplemental_unmerged_16S_reads": f"{rescued_reads:g}",
        "samples_with_supplemental_reads": sum(
            any(retained_counts[feature_id].get(sample, 0) > 0 for feature_id in rescued_ids)
            for sample in samples
        ),
        "supplemental_chloroplast_features": len(chloroplast_ids),
        "supplemental_chloroplast_reads": f"{sum(sum(retained_counts[feature_id].values()) for feature_id in chloroplast_ids):g}",
        "interpretation": "Experimental linked-pair supplement; excluded from canonical abundance and richness outputs",
    }
    return rows, summary, samples


def write_outputs(table_path, summary_path, rows, summary, samples):
    table_columns = [
        "Feature_ID",
        "Sequence_Representation",
        "Linked_Sequence",
        *AUDIT_COLUMNS[1:],
        *samples,
    ]
    for path, content, columns in (
        (table_path, rows, table_columns),
        (summary_path, [summary], list(summary)),
    ):
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(content)


def run_from_snakemake(snakemake_object):
    rows, summary, samples = build_supplement(
        read_biom_tsv(snakemake_object.input.retained_table_tsv),
        read_biom_tsv(snakemake_object.input.standard_table_tsv),
        read_fasta(snakemake_object.input.sequence_fasta),
        read_taxonomy(snakemake_object.input.silva_export),
        read_taxonomy(snakemake_object.input.pr2_export),
        float(snakemake_object.params.min_confidence),
    )
    write_outputs(
        snakemake_object.output.table,
        snakemake_object.output.summary,
        rows,
        summary,
        samples,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
