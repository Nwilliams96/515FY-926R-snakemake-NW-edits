"""Collect the forward-only sensitivity-analysis products for download."""

import shutil
from pathlib import Path


def run(snakemake_object):
    output_dir = Path(str(snakemake_object.output[0]))
    output_dir.mkdir(parents=True, exist_ok=True)

    copies = {
        "forward_only_ASV_comparison.tsv": snakemake_object.input.asv_comparison,
        "taxon_presence_comparison.tsv": snakemake_object.input.taxon_comparison,
        "forward_only_unique_ASVs.fasta": snakemake_object.input.unique_fasta,
        "forward_only_summary.tsv": snakemake_object.input.summary,
        "forward_only_feature_table.tsv": snakemake_object.input.forward_table,
        "forward_only_sequences.fasta": snakemake_object.input.forward_sequences,
        "forward_only_taxonomy.tsv": snakemake_object.input.forward_taxonomy,
        "forward_only_denoising_stats.tsv": snakemake_object.input.forward_stats,
        "PR2-plastid-routing-audit.tsv": snakemake_object.input.pr2_audit,
        "PR2-plastid-routing-summary.tsv": snakemake_object.input.pr2_summary,
        "SILVA_species_blast_hits.tsv": snakemake_object.input.blast_hits,
        "insert_length_estimates.tsv": snakemake_object.input.estimates,
        "insert_length_summary.tsv": snakemake_object.input.length_summary,
    }
    for name, source in copies.items():
        shutil.copy2(str(source), output_dir / name)

    maximum_length = (
        int(snakemake_object.params.trunc_len_f)
        + int(snakemake_object.params.trunc_len_r)
        - int(snakemake_object.params.min_overlap)
    )
    readme = f"""Forward-only 16S sensitivity test
================================

Study: {snakemake_object.params.study_name}

Purpose
-------
This analysis asks whether taxa are lost from the canonical paired-end result
because of reverse-read quality or the forward/reverse overlap requirement. It
reuses the primer-trimmed, BBsplit-assigned prokaryotic R1 FASTQs and leaves the
canonical pipeline results unchanged.

Configured geometry
-------------------
Forward truncation length: {snakemake_object.params.trunc_len_f}
Canonical reverse truncation length: {snakemake_object.params.trunc_len_r}
Canonical minimum overlap: {snakemake_object.params.min_overlap}
Maximum mergeable post-primer insert: {maximum_length}

Key files
---------
forward_only_summary.tsv
    Overall forward-only recovery and exact-prefix comparison with paired ASVs.
forward_only_ASV_comparison.tsv
    One row per forward ASV. Exact_prefix_in_paired_ASV is the strongest direct
    test of whether that forward sequence was represented in a canonical merged
    ASV. Forward_only_unique_ASV means no exact paired-ASV prefix was found.
taxon_presence_comparison.tsv
    Read and ASV totals at matching SILVA ranks in the forward and paired runs.
insert_length_estimates.tsv
    Reference-based length estimates for forward-only-unique ASVs. A false value
    in Length_mergeable_at_configured_min_overlap supports length-dependent loss.
insert_length_summary.tsv
    Counts and read abundance associated with reference-predicted long inserts.
forward_only_denoising_stats.tsv
    Per-sample DADA2 input, filtering, denoising, and non-chimeric counts.

Interpretation cautions
-----------------------
Forward-only ASVs are not interchangeable with paired-end ASVs and must not be
merged into the canonical feature table. Different full amplicons can share the
same forward sequence. Insert lengths are inferred from close SILVA species
references; they are hypotheses rather than direct measurements. An ASV without
a qualifying reference or resolvable primer sites is retained in the estimates
table with a diagnostic status instead of being assigned a length.
"""
    Path(str(snakemake_object.output.readme)).write_text(readme, encoding="utf-8")


if "snakemake" in globals():
    run(snakemake)
