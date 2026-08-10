"""Generate a self-contained HTML summary for a completed amplicon run."""

import base64
import csv
import html
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


PALETTE = [
    "#2563eb", "#0d9488", "#f59e0b", "#dc2626", "#7c3aed",
    "#0891b2", "#65a30d", "#db2777", "#4f46e5", "#ea580c",
    "#64748b",
]

TAXONOMY_RANKS = (
    "Domain", "Supergroup", "Division", "Subdivision", "Phylum", "Class",
    "Order", "Family", "Genus", "Species", "ProPortal_ASV_Ecotype",
)

TAXONOMY_EXPLORER_JS = r"""
(() => {
  const source = document.getElementById("taxonomy-explorer-data");
  if (!source) return;
  const data = JSON.parse(source.textContent);
  const fieldSelect = document.getElementById("taxonomy-metadata-field");
  const valueSelect = document.getElementById("taxonomy-metadata-value");
  const rankSelect = document.getElementById("taxonomy-rank");
  const summary = document.getElementById("taxonomy-filter-summary");
  const chart = document.getElementById("taxonomy-explorer-chart");
  const tableBody = document.getElementById("taxonomy-explorer-body");
  const numberFormat = new Intl.NumberFormat();

  function appendOption(select, value, label) {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = label;
    select.appendChild(option);
  }

  appendOption(fieldSelect, "", "All samples");
  data.metadataFields.forEach(field => appendOption(fieldSelect, field, field));
  data.ranks.forEach(rank => appendOption(
    rankSelect,
    rank,
    rank === "ProPortal_ASV_Ecotype" ? "ProPortal ASV ecotype" : rank
  ));
  if (data.ranks.includes("Domain")) rankSelect.value = "Domain";

  function refreshMetadataValues() {
    const field = fieldSelect.value;
    valueSelect.replaceChildren();
    appendOption(valueSelect, "", "All values");
    valueSelect.disabled = !field;
    if (!field) return;
    const values = new Set();
    Object.values(data.samples).forEach(sample => {
      const value = sample.metadata[field];
      if (value !== undefined && value !== null && String(value).trim() !== "") {
        values.add(String(value));
      }
    });
    Array.from(values).sort((a, b) => a.localeCompare(b, undefined, {numeric: true}))
      .forEach(value => appendOption(valueSelect, value, value));
  }

  function selectedSamples() {
    const field = fieldSelect.value;
    const value = valueSelect.value;
    return Object.entries(data.samples).filter(([, sample]) => {
      if (!field || !value) return true;
      return String(sample.metadata[field] ?? "") === value;
    });
  }

  function render() {
    const rank = rankSelect.value;
    const selected = selectedSamples();
    const totals = new Map();
    const detections = new Map();
    selected.forEach(([, sample]) => {
      const counts = sample.taxonomy[rank] || {};
      Object.entries(counts).forEach(([taxon, abundance]) => {
        const numeric = Number(abundance) || 0;
        if (numeric <= 0) return;
        totals.set(taxon, (totals.get(taxon) || 0) + numeric);
        detections.set(taxon, (detections.get(taxon) || 0) + 1);
      });
    });
    const rows = Array.from(totals.entries())
      .sort((a, b) => b[1] - a[1])
      .slice(0, 20);
    const grandTotal = Array.from(totals.values()).reduce((sum, value) => sum + value, 0);
    const maxValue = rows.length ? rows[0][1] : 1;

    const field = fieldSelect.value;
    const value = valueSelect.value;
    const filterText = field && value ? `${field} = ${value}` : "all sample metadata";
    summary.textContent = `${selected.length} sample(s) included · ${filterText} · ${rank}`;

    chart.replaceChildren();
    tableBody.replaceChildren();
    if (!rows.length) {
      const empty = document.createElement("p");
      empty.className = "small-muted";
      empty.textContent = "No positive-abundance taxonomy records match this selection.";
      chart.appendChild(empty);
      return;
    }

    rows.slice(0, 12).forEach(([taxon, abundance]) => {
      const row = document.createElement("div");
      row.className = "taxonomy-bar-row";
      const label = document.createElement("div");
      label.className = "taxonomy-bar-label";
      label.textContent = taxon;
      const track = document.createElement("div");
      track.className = "taxonomy-bar-track";
      const fill = document.createElement("div");
      fill.className = "taxonomy-bar-fill";
      fill.style.width = `${100 * abundance / maxValue}%`;
      track.appendChild(fill);
      const count = document.createElement("div");
      count.className = "taxonomy-bar-count";
      count.textContent = numberFormat.format(Math.round(abundance));
      row.append(label, track, count);
      chart.appendChild(row);
    });

    rows.forEach(([taxon, abundance], index) => {
      const tr = document.createElement("tr");
      const values = [
        index + 1,
        taxon,
        numberFormat.format(Math.round(abundance)),
        grandTotal ? `${(100 * abundance / grandTotal).toFixed(2)}%` : "—",
        detections.get(taxon) || 0,
      ];
      values.forEach(value => {
        const td = document.createElement("td");
        td.textContent = value;
        tr.appendChild(td);
      });
      tableBody.appendChild(tr);
    });
  }

  fieldSelect.addEventListener("change", () => {
    refreshMetadataValues();
    render();
  });
  valueSelect.addEventListener("change", render);
  rankSelect.addEventListener("change", render);
  refreshMetadataValues();
  render();
})();
"""


def esc(value):
    return html.escape(str(value if value is not None else ""))


def number(value):
    if value is None or value == "":
        return None
    try:
        return float(str(value).replace(",", "").strip())
    except (TypeError, ValueError):
        return None


def fmt_count(value):
    return "—" if value is None else f"{int(round(value)):,}"


def fmt_percent(value):
    return "—" if value is None else f"{100 * value:.1f}%"


def fmt_loss(before, after):
    if before is None or after is None:
        return "—"
    lost = max(0, before - after)
    fraction = lost / before if before else None
    return f"{fmt_count(lost)} ({fmt_percent(fraction)})"


def normalize_sample_id(value):
    return str(value).strip().replace("_", "-") if value is not None else ""


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def first_value(row, names):
    normalized = {re.sub(r"[^a-z0-9]", "", key.lower()): value for key, value in row.items()}
    for name in names:
        key = re.sub(r"[^a-z0-9]", "", name.lower())
        if key in normalized:
            return normalized[key]
    return None


def parse_dada2_stats(path):
    result = {}
    for row in read_tsv(path):
        sample = first_value(row, ["sample-id", "sampleid", "sample"])
        input_reads = number(first_value(row, ["input", "input reads"]))
        filtered_reads = number(first_value(row, ["filtered", "filtered reads"]))
        denoised_reads = number(first_value(row, ["denoised", "denoised reads"]))
        merged_reads = number(first_value(row, ["merged", "merged reads"]))
        final_reads = number(
            first_value(row, ["non-chimeric", "non_chimeric", "nonchimeric"])
        )
        if sample and input_reads is not None and final_reads is not None:
            result[normalize_sample_id(sample)] = {
                "input": input_reads,
                "filtered": filtered_reads,
                "denoised": denoised_reads,
                "merged": merged_reads,
                "final": final_reads,
                "retention": final_reads / input_reads if input_reads else None,
            }
    return result


def aggregate_dada2_stats(stats):
    result = {}
    for stage in ("input", "filtered", "denoised", "merged", "final"):
        values = [row.get(stage) for row in stats.values() if row.get(stage) is not None]
        result[stage] = sum(values) if values else None
    result["retention"] = (
        result["final"] / result["input"]
        if result["input"] and result["final"] is not None
        else None
    )
    return result


def flatten_config(config, prefix=""):
    rows = []
    for key in sorted(config):
        value = config[key]
        label = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            rows.extend(flatten_config(value, label))
        elif isinstance(value, (list, tuple)):
            rows.append((label, ", ".join(map(str, value))))
        else:
            rows.append((label, value))
    return rows


def effective_dada2_parameters(config):
    """Return every user-facing DADA2 setting, including legacy defaults."""
    prok = {
        "max_ee_f": 2.0,
        "max_ee_r": 2.0,
        "trunc_q": 2,
        "min_overlap": 12,
        "pooling_method": "independent",
        "chimera_method": "consensus",
        "min_fold_parent_over_abundance": 1.0,
        "n_reads_learn": 1000000,
    }
    euk = {
        "max_ee": 2.0,
        "trunc_q": 0,
        "pooling_method": "independent",
        "chimera_method": "consensus",
        "min_fold_parent_over_abundance": 1.0,
        "n_reads_learn": 1000000,
    }
    dada2 = config.get("dada2", {}) or {}
    if isinstance(dada2.get("prokaryotes"), dict):
        prok.update(dada2["prokaryotes"])
    if isinstance(dada2.get("eukaryotes"), dict):
        euk.update(dada2["eukaryotes"])
    trunc = config.get("trunclens", {}) or {}

    return [
        ("16S paired", "trunc_len_f", trunc.get("truncR1", "—"), "Forward bases retained"),
        ("16S paired", "trunc_len_r", trunc.get("truncR2", "—"), "Reverse bases retained"),
        ("16S paired", "max_ee_f", prok["max_ee_f"], "Maximum expected errors in a forward read"),
        ("16S paired", "max_ee_r", prok["max_ee_r"], "Maximum expected errors in a reverse read"),
        ("16S paired", "trunc_q", prok["trunc_q"], "Quality score that triggers read truncation"),
        ("16S paired", "min_overlap", prok["min_overlap"], "Minimum overlap required to merge a read pair"),
        ("16S paired", "pooling_method", prok["pooling_method"], "Sample pooling used during ASV inference"),
        ("16S paired", "chimera_method", prok["chimera_method"], "Chimera-detection strategy"),
        ("16S paired", "min_fold_parent_over_abundance", prok["min_fold_parent_over_abundance"], "Minimum parent abundance used for chimera detection"),
        ("16S paired", "n_reads_learn", prok["n_reads_learn"], "Reads used to train the error model"),
        ("18S concatenated", "R1_length_before_concatenation", trunc.get("truncR1", "—"), "Forward bases retained before concatenation"),
        ("18S concatenated", "R2_length_before_concatenation", trunc.get("truncR2", "—"), "Reverse bases retained before concatenation"),
        ("18S concatenated", "max_ee", euk["max_ee"], "Maximum expected errors in a concatenated read"),
        ("18S concatenated", "trunc_q", euk["trunc_q"], "Quality score that triggers read truncation"),
        ("18S concatenated", "pooling_method", euk["pooling_method"], "Sample pooling used during ASV inference"),
        ("18S concatenated", "chimera_method", euk["chimera_method"], "Chimera-detection strategy"),
        ("18S concatenated", "min_fold_parent_over_abundance", euk["min_fold_parent_over_abundance"], "Minimum parent abundance used for chimera detection"),
        ("18S concatenated", "n_reads_learn", euk["n_reads_learn"], "Reads used to train the error model"),
    ]


def svg_composition(sample_totals, title):
    samples = list(sample_totals)
    categories = []
    category_totals = Counter()
    for values in sample_totals.values():
        category_totals.update(values)
    categories = [name for name, _ in category_totals.most_common(9)]
    if len(category_totals) > len(categories):
        categories.append("Other")
    width = max(760, 80 + 42 * len(samples))
    longest_label = max((len(str(sample)) for sample in samples), default=0)
    label_space = min(145, max(80, longest_label * 6))
    height = 340 + label_space
    plot_x, plot_y, plot_w, plot_h = 65, 30, width - 100, 300
    bar_w = min(30, plot_w / max(1, len(samples)) * 0.72)
    step = plot_w / max(1, len(samples))
    pieces = [f'<div class="chart-scroll"><svg viewBox="0 0 {width} {height}" style="min-width:{width}px" role="img" aria-label="{esc(title)}">']
    for tick in range(0, 101, 25):
        y = plot_y + plot_h * (1 - tick / 100)
        pieces.append(f'<line x1="{plot_x}" y1="{y}" x2="{plot_x + plot_w}" y2="{y}" class="grid"/>')
        pieces.append(f'<text x="{plot_x - 8}" y="{y + 4}" text-anchor="end" class="svg-label">{tick}%</text>')
    for index, sample in enumerate(samples):
        values = sample_totals[sample]
        total = sum(values.values()) or 1
        x = plot_x + step * index + (step - bar_w) / 2
        y_bottom = plot_y + plot_h
        for cat_index, category in enumerate(categories):
            value = values.get(category, 0)
            if category == "Other":
                value = sum(v for k, v in values.items() if k not in categories)
            segment_h = plot_h * value / total
            y_bottom -= segment_h
            pieces.append(
                f'<rect x="{x:.1f}" y="{y_bottom:.1f}" width="{bar_w:.1f}" height="{segment_h:.1f}" '
                f'fill="{PALETTE[cat_index % len(PALETTE)]}"><title>{esc(sample)} — {esc(category)}: {100 * value / total:.1f}%</title></rect>'
            )
        pieces.append(f'<text transform="translate({x + bar_w / 2:.1f},{plot_y + plot_h + 10}) rotate(55)" class="svg-label">{esc(sample)}</text>')
    pieces.append('</svg><div class="chart-legend" aria-label="Domain legend">')
    for index, category in enumerate(categories):
        pieces.append(
            f'<span class="legend-item"><span class="legend-swatch" '
            f'style="background:{PALETTE[index % len(PALETTE)]}"></span>{esc(category)}</span>'
        )
    pieces.append('</div></div>')
    return "".join(pieces)


def svg_top_taxa(taxa_totals, limit=12):
    data = taxa_totals.most_common(limit)
    width, height = 900, max(260, 46 * len(data) + 50)
    left, right = 220, 90
    plot_w = width - left - right
    max_value = max((value for _, value in data), default=1) or 1
    pieces = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Most abundant taxa">']
    for index, (taxon, value) in enumerate(data):
        y = 24 + index * 42
        bar_w = plot_w * value / max_value
        pieces.append(f'<text x="{left - 12}" y="{y + 21}" text-anchor="end" class="svg-label">{esc(taxon)}</text>')
        pieces.append(f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="28" rx="4" fill="#2563eb"><title>{esc(taxon)}: {fmt_count(value)}</title></rect>')
        pieces.append(f'<text x="{left + bar_w + 8:.1f}" y="{y + 20}" class="svg-total">{fmt_count(value)}</text>')
    pieces.append('</svg>')
    return "".join(pieces)


def svg_assignment_counts(assignment_totals):
    order = [
        "Prokaryotic 16S",
        "Eukaryotic 18S",
        "Chloroplast 16S",
        "Mitochondrial 16S",
        "Unassigned",
        "Other",
    ]
    data = [(category, assignment_totals.get(category, 0)) for category in order]
    data = [(category, value) for category, value in data if value > 0]
    width, height = 900, max(260, 52 * len(data) + 55)
    left, right = 220, 110
    plot_w = width - left - right
    max_value = max((value for _, value in data), default=1) or 1
    pieces = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Sequence assignment counts">'
    ]
    for index, (category, value) in enumerate(data):
        y = 26 + index * 48
        bar_w = plot_w * value / max_value
        colour = PALETTE[index % len(PALETTE)]
        pieces.append(
            f'<text x="{left - 12}" y="{y + 23}" text-anchor="end" '
            f'class="svg-label">{esc(category)}</text>'
        )
        pieces.append(
            f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="32" '
            f'rx="4" fill="{colour}"><title>{esc(category)}: '
            f'{fmt_count(value)}</title></rect>'
        )
        pieces.append(
            f'<text x="{left + bar_w + 8:.1f}" y="{y + 22}" '
            f'class="svg-total">{fmt_count(value)}</text>'
        )
    pieces.append("</svg>")
    return "".join(pieces)


def embedded_image(path):
    raw = Path(path).read_bytes()
    encoded = base64.b64encode(raw).decode("ascii")
    return f'<img class="report-figure" src="data:image/png;base64,{encoded}" alt="{esc(Path(path).stem)}">'


def taxonomy_label(row):
    for key in ("Phylum", "Division", "Supergroup", "Class", "Domain"):
        value = (row.get(key) or "").strip()
        if value and value.lower() not in {"na", "nan", "unassigned"}:
            return value
    return "Unassigned"


def taxonomy_value(row, rank):
    value = str(row.get(rank) or "").strip()
    if value and value.lower() not in {"na", "nan", "none", "unassigned"}:
        return value
    return f"Unassigned at {rank.replace('_', ' ')}"


def build_taxonomy_explorer_data(sample_rows, long_rows, abundance_column):
    """Pre-aggregate taxonomy counts for an offline, interactive HTML report."""
    metadata_fields = []
    for row in sample_rows:
        for field in row:
            normalized = re.sub(r"[^a-z0-9]", "", field.lower())
            if normalized in {"sample", "sampleid"} or field in metadata_fields:
                continue
            if any(str(candidate.get(field) or "").strip() for candidate in sample_rows):
                metadata_fields.append(field)

    metadata_by_sample = {}
    sample_order = []
    for row in sample_rows:
        sample = first_value(row, ["sample", "sample-id", "sampleid"])
        sample = normalize_sample_id(sample)
        if not sample:
            continue
        sample_order.append(sample)
        metadata_by_sample[sample] = {
            field: str(row.get(field) or "").strip() for field in metadata_fields
        }

    ranks = [rank for rank in TAXONOMY_RANKS if any(rank in row for row in long_rows)]
    taxonomy_by_sample = defaultdict(lambda: defaultdict(Counter))
    for row in long_rows:
        sample = normalize_sample_id(row.get("SampleID"))
        abundance = number(row.get(abundance_column)) or 0
        if not sample or abundance <= 0:
            continue
        if sample not in sample_order:
            sample_order.append(sample)
        for rank in ranks:
            taxonomy_by_sample[sample][rank][taxonomy_value(row, rank)] += abundance

    samples = {}
    for sample in sample_order:
        samples[sample] = {
            "metadata": metadata_by_sample.get(sample, {}),
            "taxonomy": {
                rank: dict(taxonomy_by_sample[sample].get(rank, {}))
                for rank in ranks
            },
        }
    return {"metadataFields": metadata_fields, "ranks": ranks, "samples": samples}


def sequence_assignment(row):
    lineage = " ".join(
        str(row.get(key) or "")
        for key in (
            "Domain", "Supergroup", "Division", "Subdivision", "Phylum",
            "Class", "Order", "Family", "Genus", "Species",
            "ProPortal_ASV_Ecotype",
        )
    ).lower()
    if "mitochond" in lineage:
        return "Mitochondrial 16S"

    sequence_type = re.sub(
        r"[^a-z0-9]", "", str(row.get("Sequence_Type") or "").lower()
    )
    mapping = {
        "prokaryotic16s": "Prokaryotic 16S",
        "chloroplast16s": "Chloroplast 16S",
        "eukaryote18s": "Eukaryotic 18S",
        "eukaryotic18s": "Eukaryotic 18S",
        "unassigned": "Unassigned",
    }
    if sequence_type in mapping:
        return mapping[sequence_type]

    plastid = str(row.get("plastid_16S_rRNA") or "").strip().lower()
    domain = str(row.get("Domain") or "").strip().lower()
    if plastid == "yes":
        return "Chloroplast 16S"
    if domain in {"bacteria", "archaea"}:
        return "Prokaryotic 16S"
    if domain == "eukaryota":
        return "Eukaryotic 18S"
    if not domain or domain == "unassigned":
        return "Unassigned"
    return "Other"


def quality_status(retention):
    if retention is None:
        return "not available", "neutral"
    if retention >= 0.70:
        return "good", "good"
    if retention >= 0.40:
        return "review", "warn"
    return "low", "bad"


def dada2_loss_cells(stats, paired):
    merge_loss = (
        fmt_loss(stats.get("denoised"), stats.get("merged"))
        if paired
        else '<span class="pill neutral">Not applicable</span>'
    )
    pre_chimera = stats.get("merged") if paired else stats.get("denoised")
    status, status_class = quality_status(stats.get("retention"))
    return (
        f"<td>{fmt_count(stats.get('input'))}</td>"
        f"<td>{fmt_loss(stats.get('input'), stats.get('filtered'))}</td>"
        f"<td>{fmt_loss(stats.get('filtered'), stats.get('denoised'))}</td>"
        f"<td>{merge_loss}</td>"
        f"<td>{fmt_loss(pre_chimera, stats.get('final'))}</td>"
        f"<td>{fmt_count(stats.get('final'))}</td>"
        f"<td><span class='pill {status_class}'>{fmt_percent(stats.get('retention'))} · {status}</span></td>"
    )


def render_report(config, paths, output_path):
    sample_rows = read_tsv(paths["samples"])
    sample_names = [first_value(row, ["sample", "sample-id", "sampleid"]) for row in sample_rows]
    sample_names = [normalize_sample_id(sample) for sample in sample_names if sample]
    split_rows = read_tsv(paths["split_summary"])
    split = {}
    for row in split_rows:
        sample = first_value(row, ["sample", "sample-id"])
        if not sample:
            continue
        prok = number(first_value(row, ["prok_seqs_split", "prok reads"])) or 0
        euk = number(first_value(row, ["euk_seqs_split", "euk reads"])) or 0
        split[normalize_sample_id(sample)] = {"prok": prok, "euk": euk, "total": prok + euk}

    stats16 = parse_dada2_stats(paths["stats16s"])
    stats18 = parse_dada2_stats(paths["stats18s"])
    aggregate16 = aggregate_dada2_stats(stats16)
    aggregate18 = aggregate_dada2_stats(stats18)

    long_rows = read_tsv(paths["long_data"])
    abundance_column = "Corrected_Sequence_Counts"
    if long_rows and abundance_column not in long_rows[0]:
        abundance_column = "Raw_Sequence_Counts"
    taxonomy_explorer_data = build_taxonomy_explorer_data(
        sample_rows, long_rows, abundance_column
    )
    taxonomy_explorer_json = json.dumps(
        taxonomy_explorer_data, ensure_ascii=False, separators=(",", ":")
    ).replace("</", "<\\/")
    domains_by_sample = defaultdict(Counter)
    assignment_totals = Counter()
    taxa_totals = Counter()
    taxa_samples = defaultdict(set)
    asvs = set()
    for row in long_rows:
        sample = normalize_sample_id(row.get("SampleID"))
        abundance = number(row.get(abundance_column)) or 0
        if not sample or abundance <= 0:
            continue
        domain = (row.get("Domain") or "Unassigned").strip() or "Unassigned"
        assignment = sequence_assignment(row)
        taxon = taxonomy_label(row)
        domains_by_sample[sample][domain] += abundance
        assignment_totals[assignment] += abundance
        taxa_totals[taxon] += abundance
        taxa_samples[taxon].add(sample)
        if row.get("ASV_hash"):
            asvs.add(row["ASV_hash"])

    ordered_samples = [sample for sample in sample_names if sample in set(split) | set(stats16) | set(stats18)]
    for sample in sorted((set(split) | set(stats16) | set(stats18)) - set(ordered_samples)):
        ordered_samples.append(sample)

    split_total = sum(item["total"] for item in split.values())
    final16 = sum(item["final"] for item in stats16.values())
    final18 = sum(item["final"] for item in stats18.values())
    all_retentions = [item["retention"] for item in list(stats16.values()) + list(stats18.values()) if item["retention"] is not None]
    median_retention = None
    if all_retentions:
        ordered = sorted(all_retentions)
        middle = len(ordered) // 2
        median_retention = ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2

    dada2_sample_rows = []
    for sample in ordered_samples:
        s16 = stats16.get(sample, {})
        s18 = stats18.get(sample, {})
        if s16:
            dada2_sample_rows.append(
                f"<tr><td>{esc(sample)}</td><td>16S paired</td>"
                f"{dada2_loss_cells(s16, paired=True)}</tr>"
            )
        if s18:
            dada2_sample_rows.append(
                f"<tr><td>{esc(sample)}</td><td>18S concatenated</td>"
                f"{dada2_loss_cells(s18, paired=False)}</tr>"
            )

    dada2_summary_rows = (
        f"<tr><td>16S paired</td>{dada2_loss_cells(aggregate16, paired=True)}</tr>"
        f"<tr><td>18S concatenated</td>{dada2_loss_cells(aggregate18, paired=False)}</tr>"
    )

    parameter_rows = "".join(
        f"<tr><th>{esc(key)}</th><td><code>{esc(value)}</code></td></tr>"
        for key, value in flatten_config(config)
    )
    dada2_parameter_rows = "".join(
        f"<tr><td>{esc(path)}</td><th><code>{esc(parameter)}</code></th>"
        f"<td><code>{esc(value)}</code></td><td>{esc(description)}</td></tr>"
        for path, parameter, value, description in effective_dada2_parameters(config)
    )
    top_taxa_rows = "".join(
        f"<tr><td>{index}</td><td>{esc(taxon)}</td><td>{fmt_count(value)}</td><td>{len(taxa_samples[taxon])}</td></tr>"
        for index, (taxon, value) in enumerate(taxa_totals.most_common(20), 1)
    )
    assignment_order = [
        "Prokaryotic 16S",
        "Eukaryotic 18S",
        "Chloroplast 16S",
        "Mitochondrial 16S",
        "Unassigned",
        "Other",
    ]
    assignment_total = sum(assignment_totals.values())
    total_16s = sum(
        assignment_totals[category]
        for category in (
            "Prokaryotic 16S", "Chloroplast 16S", "Mitochondrial 16S"
        )
    )
    assignment_summary = [("Total 16S", total_16s)] + [
        (category, assignment_totals[category]) for category in assignment_order
    ]
    assignment_rows = "".join(
        f"<tr><td>{esc(category)}</td><td>{fmt_count(value)}</td>"
        f"<td>{fmt_percent(value / assignment_total if assignment_total else None)}</td></tr>"
        for category, value in assignment_summary
    )
    internal_paths = [path for path in paths.get("internal_standard_figures", []) if Path(path).is_file()]
    internal_section = ""
    if internal_paths:
        internal_section = f"""
        <section id="internal-standards">
          <div class="eyebrow">Optional analysis</div><h2>Internal-standard correction</h2>
          <p>These figures are produced by the internal-standard correction step and embedded in this report.</p>
          {''.join(embedded_image(path) for path in internal_paths)}
        </section>"""

    domain_chart = svg_composition({sample: domains_by_sample[sample] for sample in ordered_samples if domains_by_sample[sample]}, "Domain composition by sample")
    assignment_chart = svg_assignment_counts(assignment_totals)
    taxa_chart = svg_top_taxa(taxa_totals)
    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
    study = config.get("studyName", "Amplicon run")
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(study)} pipeline report</title>
<style>
:root{{--ink:#172033;--muted:#667085;--line:#dfe4ea;--paper:#fff;--wash:#f4f7fb;--blue:#2563eb;--teal:#0d9488}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--wash);color:var(--ink);font:15px/1.55 Inter,ui-sans-serif,system-ui,-apple-system,sans-serif}}
header{{background:linear-gradient(130deg,#12213c,#1d4ed8);color:white;padding:56px max(24px,calc((100vw - 1180px)/2)) 48px}}
header h1{{font-size:clamp(32px,5vw,56px);line-height:1.05;margin:.15em 0}} header p{{margin:0;opacity:.82}}
nav{{position:sticky;top:0;z-index:3;background:#ffffffee;backdrop-filter:blur(9px);border-bottom:1px solid var(--line);padding:12px max(24px,calc((100vw - 1180px)/2));overflow:auto;white-space:nowrap}}
nav a{{color:#344054;text-decoration:none;margin-right:22px;font-weight:650}}
main{{max-width:1180px;margin:28px auto 70px;padding:0 22px}} section{{background:var(--paper);border:1px solid var(--line);border-radius:16px;padding:28px;margin:20px 0;box-shadow:0 10px 28px #1018280a}}
h2{{font-size:25px;margin:.1em 0 .35em}} h3{{margin-top:28px}} .eyebrow{{color:var(--blue);font-size:12px;font-weight:800;letter-spacing:.12em;text-transform:uppercase}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin:24px 0}} .card{{border:1px solid var(--line);border-radius:12px;padding:18px;background:#fbfcfe}} .card strong{{display:block;font-size:28px;line-height:1.2}} .card span{{color:var(--muted);font-size:13px}}
.table-wrap{{overflow:auto;border:1px solid var(--line);border-radius:10px}} table{{width:100%;border-collapse:collapse;white-space:nowrap}} th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid var(--line)}} thead th{{background:#f8fafc;font-size:12px;text-transform:uppercase;letter-spacing:.04em}} tbody tr:last-child td,tbody tr:last-child th{{border-bottom:0}}
.pill{{display:inline-block;padding:3px 8px;border-radius:99px;font-size:12px;font-weight:750}} .pill.good{{background:#dcfce7;color:#166534}} .pill.warn{{background:#fef3c7;color:#92400e}} .pill.bad{{background:#fee2e2;color:#991b1b}} .pill.neutral{{background:#e2e8f0;color:#475569}}
.svg-label{{fill:#475467;font-size:12px}} .svg-total{{fill:#344054;font-size:12px;font-weight:700}} .grid{{stroke:#e4e7ec;stroke-width:1}} svg{{max-width:100%;height:auto}} .chart-scroll{{overflow-x:auto}}
.chart-legend{{display:flex;flex-wrap:wrap;gap:10px 22px;align-items:center;padding:12px 10px 2px;min-width:max-content}} .legend-item{{display:inline-flex;align-items:center;gap:7px;color:#475467;font-size:12px;font-weight:650}} .legend-swatch{{display:inline-block;width:12px;height:12px;border-radius:2px;flex:none}}
.explorer-controls{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:22px 0 12px}} .explorer-control label{{display:block;margin-bottom:5px;color:#344054;font-size:13px;font-weight:750}} .explorer-control select{{width:100%;padding:10px 12px;border:1px solid #cbd5e1;border-radius:8px;background:white;color:var(--ink);font:inherit}} .explorer-control select:disabled{{background:#f1f5f9;color:#94a3b8}} .small-muted{{color:var(--muted);font-size:13px}} .taxonomy-chart{{display:grid;gap:8px;margin:20px 0 26px}} .taxonomy-bar-row{{display:grid;grid-template-columns:minmax(120px,220px) minmax(180px,1fr) 90px;gap:10px;align-items:center}} .taxonomy-bar-label{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:#344054;font-size:13px}} .taxonomy-bar-track{{height:24px;background:#e8eef7;border-radius:5px;overflow:hidden}} .taxonomy-bar-fill{{height:100%;min-width:2px;background:linear-gradient(90deg,var(--blue),var(--teal));border-radius:5px}} .taxonomy-bar-count{{text-align:right;font-variant-numeric:tabular-nums;font-size:13px;font-weight:700}}
.report-figure{{display:block;max-width:100%;height:auto;border:1px solid var(--line);border-radius:10px;margin:18px auto}} code{{white-space:normal;word-break:break-word}} .note{{background:#eff6ff;border-left:4px solid var(--blue);padding:12px 15px;border-radius:6px;color:#344054}}
@media(max-width:650px){{.taxonomy-bar-row{{grid-template-columns:minmax(90px,130px) minmax(100px,1fr) 70px}}}}
@media print{{nav{{display:none}}body{{background:white}}section{{box-shadow:none;break-inside:avoid}}}}
</style></head><body>
<header><div class="eyebrow" style="color:#bfdbfe">515Y/926R amplicon workflow</div><h1>{esc(study)}</h1><p>Pipeline summary generated {esc(generated)}</p></header>
<nav><a href="#overview">Overview</a><a href="#parameters">Parameters</a><a href="#quality">DADA2</a><a href="#composition">Composition</a><a href="#taxa">Taxa</a>{'<a href="#internal-standards">Internal standards</a>' if internal_paths else ''}</nav>
<main>
<section id="overview"><div class="eyebrow">Run at a glance</div><h2>Analysis overview</h2>
<div class="cards"><div class="card"><strong>{len(sample_names):,}</strong><span>configured samples</span></div><div class="card"><strong>{fmt_count(split_total)}</strong><span>reads assigned by 16S/18S split</span></div><div class="card"><strong>{fmt_count(final16 + final18)}</strong><span>non-chimeric reads after DADA2</span></div><div class="card"><strong>{len(asvs):,}</strong><span>observed ASVs</span></div><div class="card"><strong>{fmt_percent(median_retention)}</strong><span>median DADA2 retention</span></div></div>
<p class="note">This is a rapid quality-control summary, not a substitute for inspecting unusual samples, QIIME 2 quality visualizations, or the full result tables.</p></section>
<section id="parameters"><div class="eyebrow">Reproducibility</div><h2>Parameters used</h2><h3>Effective DADA2 settings used</h3><p>This table records the values actually applied by the pipeline. For an older config without a <code>dada2</code> block, the workflow defaults are shown.</p><div class="table-wrap"><table><thead><tr><th>Path</th><th>Parameter</th><th>Value</th><th>What it controls</th></tr></thead><tbody>{dada2_parameter_rows}</tbody></table></div><h3>Complete configuration</h3><div class="table-wrap"><table><tbody>{parameter_rows}</tbody></table></div></section>
<section id="quality"><div class="eyebrow">DADA2 processing</div><h2>Where reads were lost in DADA2</h2><p>Each loss is shown as a read count and the percentage lost from the immediately preceding stage. Filtering covers DADA2 quality filtering and truncation; denoising applies the learned error model; pair merging applies only to paired 16S reads; and the final loss is chimera removal. The 18S reads were concatenated before entering single-end DADA2, so pair merging is not applicable to that path.</p><div class="table-wrap"><table><thead><tr><th>Path</th><th>DADA2 input</th><th>Filtering loss</th><th>Denoising loss</th><th>Pair-merging loss</th><th>Chimera-removal loss</th><th>Final reads</th><th>Total retention</th></tr></thead><tbody>{dada2_summary_rows}</tbody></table></div><h3>DADA2 losses by sample</h3><p>Use this table to identify whether an individual sample loses most reads during filtering, denoising, paired-read merging, or chimera removal. Values below 40% total retention are highlighted for review; these thresholds are guides rather than automatic pass/fail criteria.</p><div class="table-wrap"><table><thead><tr><th>Sample</th><th>Path</th><th>DADA2 input</th><th>Filtering loss</th><th>Denoising loss</th><th>Pair-merging loss</th><th>Chimera-removal loss</th><th>Final reads</th><th>Total retention</th></tr></thead><tbody>{''.join(dada2_sample_rows)}</tbody></table></div></section>
<section id="composition"><div class="eyebrow">Basic bar plots</div><h2>Domain composition by sample</h2><p>Bars show relative abundance from <code>{esc(abundance_column)}</code>. Hover over a segment for its value.</p>{domain_chart}<h3>Sequence assignments</h3><p>This breakdown uses the pipeline's <code>Sequence_Type</code> field and taxonomy labels. The broad 16S total includes prokaryotic, chloroplast, and mitochondrial 16S. The figure itself uses mutually exclusive categories, so each sequence count appears in only one bar.</p><div class="cards"><div class="card"><strong>{fmt_count(total_16s)}</strong><span>total 16S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Eukaryotic 18S'])}</strong><span>eukaryotic 18S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Chloroplast 16S'])}</strong><span>chloroplast 16S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Mitochondrial 16S'])}</strong><span>mitochondrial 16S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Unassigned'])}</strong><span>unassigned</span></div></div>{assignment_chart}<h3>Sequence-assignment counts</h3><div class="table-wrap"><table><thead><tr><th>Assignment</th><th>Sequence abundance</th><th>Share of all assignments</th></tr></thead><tbody>{assignment_rows}</tbody></table></div><p class="note"><strong>Total 16S</strong> is a summary row and overlaps its three 16S subcategories; the remaining rows and the figure are mutually exclusive.</p></section>
<section id="taxa"><div class="eyebrow">Taxonomic summary</div><h2>Interactive taxonomy explorer</h2><p>Use a populated field from <code>samples.tsv</code> to examine a sample group, then choose any taxonomy level available in the formatted results. Counts use <code>{esc(abundance_column)}</code>. The chart shows the 12 most abundant taxa and the table shows the top 20.</p><div class="explorer-controls"><div class="explorer-control"><label for="taxonomy-metadata-field">Sample metadata field</label><select id="taxonomy-metadata-field"></select></div><div class="explorer-control"><label for="taxonomy-metadata-value">Metadata value</label><select id="taxonomy-metadata-value"></select></div><div class="explorer-control"><label for="taxonomy-rank">Taxonomy level</label><select id="taxonomy-rank"></select></div></div><p id="taxonomy-filter-summary" class="small-muted" aria-live="polite"></p><div id="taxonomy-explorer-chart" class="taxonomy-chart" aria-label="Filtered taxonomic abundance chart"></div><h3>Top taxa table</h3><div class="table-wrap"><table><thead><tr><th>#</th><th>Taxon</th><th>Total abundance</th><th>Relative abundance</th><th>Samples detected</th></tr></thead><tbody id="taxonomy-explorer-body"></tbody></table></div><noscript><p class="note">Interactive filters require JavaScript. This static summary uses all samples and the first informative SILVA or PR2 rank.</p>{taxa_chart}<div class="table-wrap"><table><thead><tr><th>#</th><th>Taxon</th><th>Total abundance</th><th>Samples detected</th></tr></thead><tbody>{top_taxa_rows}</tbody></table></div></noscript></section>
{internal_section}
</main><script type="application/json" id="taxonomy-explorer-data">{taxonomy_explorer_json}</script><script>{TAXONOMY_EXPLORER_JS}</script></body></html>"""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")


def run_from_snakemake(snakemake_object):
    paths = {
        "samples": str(snakemake_object.input.samples),
        "split_summary": str(snakemake_object.input.split_summary),
        "stats16s": str(snakemake_object.input.stats16s),
        "stats18s": str(snakemake_object.input.stats18s),
        "long_data": str(snakemake_object.input.long_data),
        "internal_standard_figures": list(
            getattr(snakemake_object.input, "internal_standard_figures", []) or []
        ),
    }
    render_report(dict(snakemake_object.config), paths, str(snakemake_object.output.html))


if "snakemake" in globals():
    run_from_snakemake(snakemake)
