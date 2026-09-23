"""Compare forward-only ASVs with the canonical paired-end 16S result."""

import csv
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


RANKS = {
    "d__": "Domain",
    "k__": "Kingdom",
    "p__": "Phylum",
    "c__": "Class",
    "o__": "Order",
    "f__": "Family",
    "g__": "Genus",
    "s__": "Species",
}


def run_command(parts, log_handle):
    subprocess.run(
        [str(part) for part in parts],
        check=True,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
    )


def export_artifact(artifact, output_dir, log_handle):
    run_command(
        ["qiime", "tools", "export", "--input-path", artifact,
         "--output-path", output_dir],
        log_handle,
    )


def read_fasta(path):
    records = {}
    identifier = None
    sequence = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence).upper()
                identifier = line[1:].split()[0]
                sequence = []
            elif identifier is not None:
                sequence.append(line)
    if identifier is not None:
        records[identifier] = "".join(sequence).upper()
    return records


def read_table(path):
    counts = {}
    sample_ids = []
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#OTU ID") or row[0] == "Feature ID":
                sample_ids = row[1:]
                continue
            if row[0].startswith("#"):
                continue
            values = [float(value or 0) for value in row[1:]]
            counts[row[0]] = dict(zip(sample_ids, values))
    return counts, sample_ids


def read_taxonomy(path):
    result = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature_id = (row.get("Feature ID") or row.get("FeatureID") or "").strip()
            if feature_id:
                result[feature_id] = {
                    "taxonomy": (row.get("Taxon") or row.get("Taxonomy") or "").strip(),
                    "confidence": (row.get("Confidence") or "").strip(),
                }
    return result


def named_ranks(taxonomy):
    found = {}
    for token in taxonomy.split(";"):
        token = token.strip()
        prefix = token[:3]
        value = token[3:].strip() if prefix in RANKS else ""
        if value:
            found[RANKS[prefix]] = value
    return found


def write_tsv(path, columns, rows):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, identifiers, sequences):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for feature_id in identifiers:
            handle.write(f">{feature_id}\n{sequences[feature_id]}\n")


def aggregate_taxa(counts, taxonomy):
    result = defaultdict(lambda: {"reads": 0.0, "features": set()})
    for feature_id, sample_counts in counts.items():
        total = sum(sample_counts.values())
        for rank, taxon in named_ranks(
            taxonomy.get(feature_id, {}).get("taxonomy", "")
        ).items():
            key = (rank, taxon)
            result[key]["reads"] += total
            result[key]["features"].add(feature_id)
    return result


def format_number(value):
    return str(int(value)) if float(value).is_integer() else f"{value:.10g}"


def run(snakemake_object):
    log_path = Path(str(snakemake_object.log[0]))
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="forward-only-16S-") as temp_dir:
        root = Path(temp_dir)
        forward = root / "forward"
        paired = root / "paired"
        for directory in (
            forward / "table", forward / "sequences", forward / "taxonomy",
            forward / "silva-taxonomy", forward / "stats", paired / "table",
            paired / "sequences", paired / "taxonomy", paired / "silva-taxonomy",
        ):
            directory.mkdir(parents=True, exist_ok=True)

        with log_path.open("w", encoding="utf-8") as log_handle:
            export_artifact(snakemake_object.input.forward_table, forward / "table", log_handle)
            export_artifact(snakemake_object.input.forward_sequences, forward / "sequences", log_handle)
            export_artifact(snakemake_object.input.forward_taxonomy, forward / "taxonomy", log_handle)
            export_artifact(snakemake_object.input.forward_silva_taxonomy, forward / "silva-taxonomy", log_handle)
            export_artifact(snakemake_object.input.forward_stats, forward / "stats", log_handle)
            export_artifact(snakemake_object.input.paired_table, paired / "table", log_handle)
            export_artifact(snakemake_object.input.paired_sequences, paired / "sequences", log_handle)
            export_artifact(snakemake_object.input.paired_taxonomy, paired / "taxonomy", log_handle)
            export_artifact(snakemake_object.input.paired_silva_taxonomy, paired / "silva-taxonomy", log_handle)

            forward_table_tsv = forward / "feature-table.tsv"
            paired_table_tsv = paired / "feature-table.tsv"
            run_command(
                ["biom", "convert", "-i", forward / "table/feature-table.biom",
                 "-o", forward_table_tsv, "--to-tsv"],
                log_handle,
            )
            run_command(
                ["biom", "convert", "-i", paired / "table/feature-table.biom",
                 "-o", paired_table_tsv, "--to-tsv"],
                log_handle,
            )

        forward_counts, forward_samples = read_table(forward_table_tsv)
        paired_counts, _ = read_table(paired_table_tsv)
        forward_sequences = read_fasta(forward / "sequences/dna-sequences.fasta")
        paired_sequences = read_fasta(paired / "sequences/dna-sequences.fasta")
        forward_taxonomy = read_taxonomy(forward / "taxonomy/taxonomy.tsv")
        paired_taxonomy = read_taxonomy(paired / "taxonomy/taxonomy.tsv")
        forward_silva_taxonomy = read_taxonomy(forward / "silva-taxonomy/taxonomy.tsv")
        paired_silva_taxonomy = read_taxonomy(paired / "silva-taxonomy/taxonomy.tsv")

        prefix_indexes = {}
        for length in {len(seq) for seq in forward_sequences.values()}:
            index = defaultdict(list)
            for paired_id, paired_seq in paired_sequences.items():
                if len(paired_seq) >= length:
                    index[paired_seq[:length]].append(paired_id)
            prefix_indexes[length] = index

        paired_taxonomy_values = {
            record["taxonomy"] for record in paired_taxonomy.values()
            if record["taxonomy"]
        }
        comparison_rows = []
        unique_ids = []
        unique_reads = 0.0
        all_reads = 0.0
        for feature_id in sorted(
            forward_sequences,
            key=lambda value: (-sum(forward_counts.get(value, {}).values()), value),
        ):
            sequence = forward_sequences[feature_id]
            matches = prefix_indexes[len(sequence)].get(sequence, [])
            sample_counts = forward_counts.get(feature_id, {})
            total = sum(sample_counts.values())
            all_reads += total
            taxonomy_record = forward_taxonomy.get(
                feature_id, {"taxonomy": "", "confidence": ""}
            )
            is_unique = not matches
            if is_unique:
                unique_ids.append(feature_id)
                unique_reads += total
            comparison_rows.append(
                {
                    "Forward_ASV_ID": feature_id,
                    "Forward_sequence_length": len(sequence),
                    "Total_forward_reads": format_number(total),
                    "Samples_detected": sum(value > 0 for value in sample_counts.values()),
                    "Exact_prefix_in_paired_ASV": str(bool(matches)).lower(),
                    "Matching_paired_ASV_count": len(matches),
                    "Matching_paired_ASV_IDs": ";".join(sorted(matches)),
                    "Forward_only_unique_ASV": str(is_unique).lower(),
                    "Taxonomy": taxonomy_record["taxonomy"],
                    "Taxonomy_confidence": taxonomy_record["confidence"],
                    "Same_taxonomy_present_in_paired": str(
                        bool(taxonomy_record["taxonomy"])
                        and taxonomy_record["taxonomy"] in paired_taxonomy_values
                    ).lower(),
                }
            )

        comparison_columns = list(comparison_rows[0]) if comparison_rows else [
            "Forward_ASV_ID", "Forward_sequence_length", "Total_forward_reads",
            "Samples_detected", "Exact_prefix_in_paired_ASV",
            "Matching_paired_ASV_count", "Matching_paired_ASV_IDs",
            "Forward_only_unique_ASV", "Taxonomy", "Taxonomy_confidence",
            "Same_taxonomy_present_in_paired",
        ]
        write_tsv(snakemake_object.output.asv_comparison, comparison_columns, comparison_rows)
        write_fasta(snakemake_object.output.unique_fasta, unique_ids, forward_sequences)

        forward_taxa = aggregate_taxa(forward_counts, forward_silva_taxonomy)
        paired_taxa = aggregate_taxa(paired_counts, paired_silva_taxonomy)
        taxon_rows = []
        rank_order = {rank: index for index, rank in enumerate(RANKS.values())}
        for rank, taxon in sorted(
            set(forward_taxa) | set(paired_taxa),
            key=lambda item: (rank_order.get(item[0], 99), item[1]),
        ):
            f_record = forward_taxa.get((rank, taxon), {"reads": 0.0, "features": set()})
            p_record = paired_taxa.get((rank, taxon), {"reads": 0.0, "features": set()})
            if f_record["features"] and p_record["features"]:
                status = "both"
            elif f_record["features"]:
                status = "forward_only"
            else:
                status = "paired_only"
            taxon_rows.append(
                {
                    "Rank": rank,
                    "Taxon": taxon,
                    "Presence": status,
                    "Forward_ASVs": len(f_record["features"]),
                    "Forward_reads": format_number(f_record["reads"]),
                    "Paired_ASVs": len(p_record["features"]),
                    "Paired_reads": format_number(p_record["reads"]),
                }
            )
        write_tsv(
            snakemake_object.output.taxon_comparison,
            ["Rank", "Taxon", "Presence", "Forward_ASVs", "Forward_reads",
             "Paired_ASVs", "Paired_reads"],
            taxon_rows,
        )

        summary_rows = [
            {"Metric": "forward_only_samples", "Value": len(forward_samples)},
            {"Metric": "forward_only_ASVs", "Value": len(forward_sequences)},
            {"Metric": "forward_only_reads", "Value": format_number(all_reads)},
            {"Metric": "ASVs_exactly_represented_as_paired_prefix", "Value": len(forward_sequences) - len(unique_ids)},
            {"Metric": "ASVs_not_exactly_represented_as_paired_prefix", "Value": len(unique_ids)},
            {"Metric": "reads_in_ASVs_not_exactly_represented_as_paired_prefix", "Value": format_number(unique_reads)},
            {"Metric": "paired_end_ASVs", "Value": len(paired_sequences)},
            {"Metric": "paired_end_reads", "Value": format_number(sum(sum(row.values()) for row in paired_counts.values()))},
        ]
        write_tsv(snakemake_object.output.summary, ["Metric", "Value"], summary_rows)

        shutil.copyfile(forward_table_tsv, snakemake_object.output.forward_table)
        shutil.copyfile(
            forward / "sequences/dna-sequences.fasta",
            snakemake_object.output.forward_sequences,
        )
        shutil.copyfile(
            forward / "taxonomy/taxonomy.tsv",
            snakemake_object.output.forward_taxonomy,
        )
        shutil.copyfile(
            forward / "stats/stats.tsv",
            snakemake_object.output.forward_stats,
        )


if "snakemake" in globals():
    run(snakemake)
