"""Conservatively replace SILVA taxonomy with confident PR2 plastid calls."""

import csv
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


AUDIT_COLUMNS = [
    "ASV_hash",
    "SILVA_taxonomy",
    "SILVA_confidence",
    "SILVA_chloroplast",
    "PR2_taxonomy",
    "PR2_confidence",
    "PR2_plastid",
    "PR2_override_accepted",
    "Final_taxonomy",
    "Final_confidence",
    "Chloroplast_detection_source",
]


def _first(row, names):
    for name in names:
        if name in row:
            return (row.get(name) or "").strip()
    return ""


def _confidence(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_taxonomy(path):
    """Read an exported QIIME taxonomy TSV without assuming one header style."""
    rows = []
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature_id = _first(row, ["Feature ID", "FeatureID", "feature-id", "ID"])
            if not feature_id or feature_id.startswith("#"):
                continue
            rows.append(
                {
                    "feature_id": feature_id,
                    "taxonomy": _first(row, ["Taxon", "Taxonomy", "taxon"]),
                    "confidence_text": _first(row, ["Confidence", "confidence"]),
                    "confidence": _confidence(
                        _first(row, ["Confidence", "confidence"])
                    ),
                }
            )
    return rows


def is_silva_chloroplast(taxonomy):
    return bool(
        re.search(r"(?:^|;\s*)o__Chloroplast(?:\s*;|\s*$)", taxonomy or "", re.I)
    )


def is_pr2_plastid(taxonomy):
    return ":plas" in (taxonomy or "").lower()


def resolve_taxonomies(silva_rows, pr2_rows, min_confidence=0.7):
    """Return final QIIME taxonomy rows, a per-ASV audit, and summary counts."""
    pr2_by_id = {row["feature_id"]: row for row in pr2_rows}
    final_rows = []
    audit_rows = []
    counts = Counter()

    for silva in silva_rows:
        feature_id = silva["feature_id"]
        pr2 = pr2_by_id.get(
            feature_id,
            {
                "taxonomy": "",
                "confidence_text": "",
                "confidence": None,
            },
        )
        silva_chloroplast = is_silva_chloroplast(silva["taxonomy"])
        pr2_plastid = is_pr2_plastid(pr2["taxonomy"])
        pr2_confidence = pr2["confidence"]
        accept_pr2 = bool(
            pr2_plastid
            and pr2_confidence is not None
            and pr2_confidence >= min_confidence
        )

        if accept_pr2:
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

        final_rows.append(
            {
                "Feature ID": feature_id,
                "Taxon": final_taxonomy,
                "Confidence": final_confidence,
            }
        )
        audit_rows.append(
            {
                "ASV_hash": feature_id,
                "SILVA_taxonomy": silva["taxonomy"],
                "SILVA_confidence": silva["confidence_text"],
                "SILVA_chloroplast": str(silva_chloroplast).lower(),
                "PR2_taxonomy": pr2["taxonomy"],
                "PR2_confidence": pr2["confidence_text"],
                "PR2_plastid": str(pr2_plastid).lower(),
                "PR2_override_accepted": str(accept_pr2).lower(),
                "Final_taxonomy": final_taxonomy,
                "Final_confidence": final_confidence,
                "Chloroplast_detection_source": source,
            }
        )
        counts[source] += 1
        counts["SILVA_chloroplast"] += int(silva_chloroplast)
        counts["PR2_plastid_above_threshold"] += int(accept_pr2)
        counts["PR2_plastid_any_confidence"] += int(pr2_plastid)

    summary = {
        "PR2_min_confidence": f"{min_confidence:g}",
        "total_16S_ASVs": len(silva_rows),
        "SILVA_chloroplast_ASVs": counts["SILVA_chloroplast"],
        "PR2_plastid_ASVs_any_confidence": counts["PR2_plastid_any_confidence"],
        "PR2_plastid_ASVs_accepted": counts["PR2_plastid_above_threshold"],
        "PR2_confirmed_SILVA_ASVs": counts["PR2_confirmed_SILVA"],
        "PR2_rescued_ASVs": counts["PR2_rescue"],
        "SILVA_chloroplast_unconfirmed_ASVs": counts[
            "SILVA_chloroplast_unconfirmed"
        ],
        "PR2_low_confidence_not_accepted_ASVs": counts[
            "PR2_low_confidence_not_accepted"
        ],
        "SILVA_retained_ASVs": counts["SILVA_retained"],
    }
    return final_rows, audit_rows, summary


def write_tsv(path, rows, columns):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _run(command, log_handle):
    subprocess.run(
        [str(part) for part in command],
        check=True,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
    )


def run_from_snakemake(snakemake_object):
    log_path = Path(str(snakemake_object.log[0]))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    min_confidence = float(snakemake_object.params.min_confidence)

    with tempfile.TemporaryDirectory(prefix="resolve-pr2-plastids-") as temp_dir:
        temp_root = Path(temp_dir)
        silva_export = temp_root / "silva"
        pr2_export = temp_root / "pr2"
        final_tsv = temp_root / "resolved-taxonomy.tsv"

        with log_path.open("w", encoding="utf-8") as log_handle:
            _run(
                [
                    "qiime",
                    "tools",
                    "export",
                    "--input-path",
                    snakemake_object.input.silva_taxonomy,
                    "--output-path",
                    silva_export,
                ],
                log_handle,
            )
            _run(
                [
                    "qiime",
                    "tools",
                    "export",
                    "--input-path",
                    snakemake_object.input.pr2_taxonomy,
                    "--output-path",
                    pr2_export,
                ],
                log_handle,
            )

            final_rows, audit_rows, summary = resolve_taxonomies(
                read_taxonomy(silva_export / "taxonomy.tsv"),
                read_taxonomy(pr2_export / "taxonomy.tsv"),
                min_confidence,
            )
            write_tsv(final_tsv, final_rows, ["Feature ID", "Taxon", "Confidence"])
            write_tsv(snakemake_object.output.audit, audit_rows, AUDIT_COLUMNS)
            write_tsv(snakemake_object.output.summary, [summary], list(summary))

            _run(
                [
                    "qiime",
                    "tools",
                    "import",
                    "--type",
                    "FeatureData[Taxonomy]",
                    "--input-path",
                    final_tsv,
                    "--output-path",
                    snakemake_object.output.taxonomy,
                ],
                log_handle,
            )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
