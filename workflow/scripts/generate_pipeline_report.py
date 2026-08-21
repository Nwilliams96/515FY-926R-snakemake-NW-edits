"""Generate a self-contained HTML summary for a completed amplicon run."""

import base64
import csv
import html
import json
import math
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

TAXONOMY_PLOTTING_FIELDS = (
    ("SampleID", ("sample", "sample-id", "sampleid")),
    ("Condition", ("condition", "sample condition", "sample_type", "sample type")),
    ("Latitude", ("Latitude [degrees_north]", "latitude", "lat")),
    ("Longitude", ("Longitude [degrees_east]", "longitude", "lon", "long")),
    ("Depth", ("Depth (m)", "depth", "depth_m", "depth m")),
)

TAXONOMY_EXPLORER_JS = r"""
(() => {
  const source = document.getElementById("taxonomy-explorer-data");
  if (!source) return;
  const data = JSON.parse(source.textContent);
  const fieldSelect = document.getElementById("taxonomy-plot-field");
  const valueSelect = document.getElementById("taxonomy-metadata-value");
  const rankSelect = document.getElementById("taxonomy-rank");
  const summary = document.getElementById("taxonomy-filter-summary");
  const chart = document.getElementById("taxonomy-explorer-chart");
  const tableBody = document.getElementById("taxonomy-explorer-body");
  const legend = document.getElementById("taxonomy-explorer-legend");
  const numberFormat = new Intl.NumberFormat();
  const palette = [
    "#2563eb", "#0d9488", "#f59e0b", "#dc2626", "#7c3aed", "#0891b2",
    "#65a30d", "#db2777", "#4f46e5", "#ea580c", "#14b8a6", "#9333ea",
    "#84cc16", "#e11d48", "#0284c7", "#d97706", "#16a34a", "#6366f1",
    "#be123c", "#475569"
  ];

  function appendOption(select, value, label) {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = label;
    select.appendChild(option);
  }

  data.metadataFields.forEach(field => appendOption(fieldSelect, field, field));
  data.ranks.forEach(rank => appendOption(
    rankSelect,
    rank,
    rank === "ProPortal_ASV_Ecotype" ? "ProPortal ASV ecotype" : rank
  ));
  if (data.ranks.includes("Phylum")) rankSelect.value = "Phylum";
  else if (data.ranks.includes("Division")) rankSelect.value = "Division";
  else if (data.ranks.includes("Domain")) rankSelect.value = "Domain";

  function refreshMetadataValues() {
    const field = fieldSelect.value;
    valueSelect.replaceChildren();
    appendOption(valueSelect, "", "All values");
    valueSelect.disabled = field === "SampleID";
    if (field === "SampleID") return;
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
    const selected = Object.entries(data.samples).filter(([, sample]) => {
      if (field === "SampleID" || !value) return true;
      return String(sample.metadata[field] ?? "") === value;
    });
    const collator = new Intl.Collator(undefined, {numeric: true, sensitivity: "base"});
    return selected.sort((a, b) => {
      const aValue = field === "SampleID" ? a[0] : String(a[1].metadata[field] ?? "");
      const bValue = field === "SampleID" ? b[0] : String(b[1].metadata[field] ?? "");
      if (!aValue && bValue) return 1;
      if (aValue && !bValue) return -1;
      const comparison = collator.compare(aValue, bValue);
      return comparison || collator.compare(a[0], b[0]);
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
    const rows = Array.from(totals.entries()).sort((a, b) => b[1] - a[1]);
    const grandTotal = Array.from(totals.values()).reduce((sum, value) => sum + value, 0);
    const displayedTaxa = rows.slice(0, 19).map(([taxon]) => taxon);
    const displayedSet = new Set(displayedTaxa);
    const hasOther = rows.some(([taxon]) => !displayedSet.has(taxon));
    const legendTaxa = hasOther ? [...displayedTaxa, "Other"] : displayedTaxa;

    const field = fieldSelect.value;
    const value = valueSelect.value;
    const filterText = field !== "SampleID" && value ? `${field} = ${value}` : "all samples";
    summary.textContent = `${selected.length} sample(s) included · ordered by ${field} · ${filterText} · ${rank}`;

    chart.replaceChildren();
    legend.replaceChildren();
    tableBody.replaceChildren();
    if (!rows.length) {
      const empty = document.createElement("p");
      empty.className = "small-muted";
      empty.textContent = "No positive-abundance taxonomy records match this selection.";
      chart.appendChild(empty);
      return;
    }

    const plot = document.createElement("div");
    plot.className = "taxonomy-stacked-plot";
    selected.forEach(([sampleId, sample]) => {
      const counts = sample.taxonomy[rank] || {};
      const total = Object.values(counts).reduce((sum, abundance) => sum + (Number(abundance) || 0), 0);
      const column = document.createElement("div");
      column.className = "taxonomy-sample-column";
      const bar = document.createElement("div");
      bar.className = "taxonomy-stacked-bar";
      displayedTaxa.forEach((taxon, index) => {
        const abundance = Number(counts[taxon]) || 0;
        if (abundance <= 0 || total <= 0) return;
        const segment = document.createElement("div");
        segment.className = "taxonomy-segment";
        segment.style.height = `${100 * abundance / total}%`;
        segment.style.background = palette[index % palette.length];
        segment.title = `${sampleId} — ${taxon}: ${(100 * abundance / total).toFixed(2)}% (${numberFormat.format(Math.round(abundance))})`;
        bar.appendChild(segment);
      });
      if (hasOther && total > 0) {
        const other = Object.entries(counts)
          .filter(([taxon]) => !displayedSet.has(taxon))
          .reduce((sum, [, abundance]) => sum + (Number(abundance) || 0), 0);
        if (other > 0) {
          const segment = document.createElement("div");
          segment.className = "taxonomy-segment";
          segment.style.height = `${100 * other / total}%`;
          segment.style.background = "#94a3b8";
          segment.title = `${sampleId} — Other: ${(100 * other / total).toFixed(2)}% (${numberFormat.format(Math.round(other))})`;
          bar.appendChild(segment);
        }
      }
      const label = document.createElement("div");
      label.className = "taxonomy-sample-label";
      const plotValue = field === "SampleID" ? sampleId : String(sample.metadata[field] ?? "");
      label.textContent = plotValue || sampleId;
      label.title = field === "SampleID" ? sampleId : `${sampleId} — ${field}: ${plotValue || "blank"}`;
      column.append(bar, label);
      plot.appendChild(column);
    });
    chart.appendChild(plot);

    legendTaxa.forEach((taxon, index) => {
      const item = document.createElement("span");
      item.className = "legend-item";
      const swatch = document.createElement("span");
      swatch.className = "legend-swatch";
      swatch.style.background = taxon === "Other" ? "#94a3b8" : palette[index % palette.length];
      item.append(swatch, document.createTextNode(taxon));
      legend.appendChild(item);
    });

    rows.slice(0, 20).forEach(([taxon, abundance], index) => {
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
        parsed = float(str(value).replace(",", "").strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def fmt_count(value):
    parsed = number(value)
    return "—" if parsed is None else f"{int(round(parsed)):,}"


def fmt_percent(value):
    parsed = number(value)
    return "—" if parsed is None else f"{100 * parsed:.1f}%"


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


def parse_cutadapt_qc(paths):
    """Read raw pair counts from the per-sample Cutadapt summary files."""
    result = {}
    count_pattern = re.compile(
        r"^\s*Total (?:read pairs|reads) processed:\s*([0-9,]+)", re.MULTILINE
    )
    for path in paths:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        match = count_pattern.search(text)
        if not match:
            continue
        sample = re.sub(r"\.qc\.txt$", "", Path(path).name)
        result[normalize_sample_id(sample)] = float(match.group(1).replace(",", ""))
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


def svg_sample_read_counts(counts, sample_order, label):
    data = [(sample, counts[sample]) for sample in sample_order if sample in counts]
    if not data:
        return '<p class="small-muted">Read counts were not available.</p>'
    width = max(760, 110 + len(data) * 54)
    height = 420
    left, right, top, bottom = 70, 25, 35, 100
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_value = max(value for _, value in data) or 1
    bar_w = min(34, plot_w / max(1, len(data)) * 0.65)
    pieces = [
        f'<div class="chart-scroll"><svg viewBox="0 0 {width} {height}" '
        f'style="min-width:{width}px" role="img" aria-label="{esc(label)}">'
    ]
    for fraction in (0, 0.25, 0.5, 0.75, 1):
        y = top + plot_h * (1 - fraction)
        pieces.append(f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" class="grid"/>')
        pieces.append(
            f'<text x="{left-8}" y="{y+4:.1f}" text-anchor="end" class="svg-label">'
            f'{fmt_count(max_value * fraction)}</text>'
        )
    slot_w = plot_w / len(data)
    for index, (sample, value) in enumerate(data):
        x = left + slot_w * (index + 0.5) - bar_w / 2
        bar_h = plot_h * value / max_value
        y = top + plot_h - bar_h
        pieces.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" '
            f'rx="3" fill="#2563eb"><title>{esc(sample)}: {fmt_count(value)}</title></rect>'
        )
        pieces.append(
            f'<text transform="translate({x + bar_w / 2:.1f},{height-bottom+12}) rotate(55)" '
            f'class="svg-label">{esc(sample)}</text>'
        )
    pieces.append("</svg></div>")
    return "".join(pieces)


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
    return (
        '<div class="figure-scroll">'
        f'<img class="report-figure" src="data:image/png;base64,{encoded}" '
        f'alt="{esc(Path(path).stem)}">'
        '</div>'
    )


def internal_method_label(column, standard_ids):
    label = re.sub(r"^Copies_", "", column)
    for index, standard_id in reversed(list(enumerate(standard_ids, 1))):
        label = label.replace(f"isd_{index}", standard_id)
    return label.replace("_mean_recovery_ratio", " mean").replace(
        "_recovery_ratio", ""
    ).replace("_", " ")


def svg_log_series(sample_order, values, series_labels, title):
    positives = [
        parsed for sample in sample_order for value in values.get(sample, {}).values()
        if (parsed := number(value)) is not None and parsed > 0
    ]
    if not positives:
        return '<p class="small-muted">No positive values were available.</p>'
    width = max(820, 150 + len(sample_order) * 48)
    height = 470
    left, right, top, bottom = 90, 30, 35, 115
    plot_w, plot_h = width - left - right, height - top - bottom
    low, high = math.log10(min(positives)), math.log10(max(positives))
    if high - low < 0.5:
        low -= 0.25
        high += 0.25
    colours = {
        label: PALETTE[index % len(PALETTE)]
        for index, label in enumerate(series_labels)
    }
    pieces = [
        f'<h4>{esc(title)}</h4><div class="chart-scroll"><svg viewBox="0 0 {width} {height}" '
        f'style="min-width:{width}px" role="img" aria-label="{esc(title)}">'
    ]
    for index in range(5):
        fraction = index / 4
        y = top + plot_h * (1 - fraction)
        tick = 10 ** (low + (high - low) * fraction)
        pieces.append(f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" class="grid"/>')
        pieces.append(f'<text x="{left-8}" y="{y+4:.1f}" text-anchor="end" class="svg-label">{fmt_count(tick)}</text>')
    x_step = plot_w / max(1, len(sample_order))
    for series_index, label in enumerate(series_labels):
        points = []
        for sample_index, sample in enumerate(sample_order):
            value = number(values.get(sample, {}).get(label))
            if value is None or value <= 0:
                continue
            x = left + x_step * (sample_index + 0.5)
            y = top + plot_h * (1 - (math.log10(value) - low) / (high - low))
            points.append((x, y, sample, value))
        if len(points) > 1:
            path = " ".join(
                ("M" if index == 0 else "L") + f" {x:.1f} {y:.1f}"
                for index, (x, y, _, _) in enumerate(points)
            )
            pieces.append(f'<path d="{path}" fill="none" stroke="{colours[label]}" stroke-width="1.5" opacity=".72"/>')
        for x, y, sample, value in points:
            pieces.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colours[label]}">'
                f'<title>{esc(sample)} — {esc(label)}: {value:.6g}</title></circle>'
            )
    for sample_index, sample in enumerate(sample_order):
        x = left + x_step * (sample_index + 0.5)
        pieces.append(
            f'<text transform="translate({x:.1f},{height-bottom+12}) rotate(55)" '
            f'class="svg-label">{esc(sample)}</text>'
        )
    pieces.append('</svg></div><div class="chart-legend">')
    for label in series_labels:
        pieces.append(
            f'<span class="legend-item"><span class="legend-swatch" '
            f'style="background:{colours[label]}"></span>{esc(label)}</span>'
        )
    pieces.append("</div>")
    return "".join(pieces)


def internal_standard_html(config, corrected_rows):
    if not corrected_rows:
        return ""
    configured_standards = config.get("intstds", [])
    if isinstance(configured_standards, dict):
        configured_standards = configured_standards.values()
    standard_ids = [str(value) for value in configured_standards]
    samples = sorted({normalize_sample_id(row.get("SampleID")) for row in corrected_rows if row.get("SampleID")})
    recovery_columns = [
        column for column in corrected_rows[0]
        if column.startswith("isd_") and column.endswith("recovery_ratio")
        or column in {"recovery_mean", "recovery_median"}
    ]
    recovery_labels = [internal_method_label(column, standard_ids) for column in recovery_columns]
    recovery_values = defaultdict(dict)
    for row in corrected_rows:
        sample = normalize_sample_id(row.get("SampleID"))
        for column, label in zip(recovery_columns, recovery_labels):
            recovery_values[sample][label] = number(row.get(column))

    copy_columns = [column for column in corrected_rows[0] if column.startswith("Copies_")]
    domain_totals = {
        column: defaultdict(lambda: defaultdict(float)) for column in copy_columns
    }
    for row in corrected_rows:
        sample = normalize_sample_id(row.get("SampleID"))
        domain = (row.get("Domain") or "Unassigned").strip() or "Unassigned"
        if domain == "Unassigned":
            continue
        for column in copy_columns:
            domain_totals[column][sample][domain] += number(row.get(column)) or 0

    recovery_chart = svg_log_series(
        samples, recovery_values, recovery_labels, "Recovery ratios by sample (log10 scale)"
    )
    domain_charts = []
    for column in copy_columns:
        method = internal_method_label(column, standard_ids)
        domains = sorted({domain for sample in domain_totals[column].values() for domain in sample})
        chart = svg_log_series(
            samples, domain_totals[column], domains,
            f"Domain copies per unit — {method} (log10 scale)",
        )
        domain_charts.append(f'<details><summary>{esc(method)}</summary>{chart}</details>')
    return recovery_chart + "<h3>Domain abundance by correction method</h3>" + "".join(domain_charts)


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
    headers = list(sample_rows[0]) if sample_rows else []
    normalized_headers = {
        re.sub(r"[^a-z0-9]", "", header.lower()): header for header in headers
    }
    metadata_sources = [("SampleID", None)]
    consumed_headers = set()
    for label, aliases in TAXONOMY_PLOTTING_FIELDS[1:]:
        source = next(
            (
                normalized_headers[re.sub(r"[^a-z0-9]", "", alias.lower())]
                for alias in aliases
                if re.sub(r"[^a-z0-9]", "", alias.lower()) in normalized_headers
            ),
            None,
        )
        if source:
            metadata_sources.append((label, source))
            consumed_headers.add(source)

    sample_headers = {
        normalized_headers.get(re.sub(r"[^a-z0-9]", "", alias.lower()))
        for alias in TAXONOMY_PLOTTING_FIELDS[0][1]
    }
    technical_headers = {
        "internal_std_normalization_factor", "units"
    }
    for header in headers:
        if (
            header not in consumed_headers
            and header not in sample_headers
            and header not in technical_headers
            and not header.endswith("_ng")
        ):
            metadata_sources.append((header, header))
    metadata_fields = [label for label, _ in metadata_sources]

    metadata_by_sample = {}
    sample_order = []
    for row in sample_rows:
        sample = first_value(row, ["sample", "sample-id", "sampleid"])
        sample = normalize_sample_id(sample)
        if not sample:
            continue
        sample_order.append(sample)
        metadata = {"SampleID": sample}
        for label, source in metadata_sources[1:]:
            metadata[label] = str(row.get(source) or "").strip()
        metadata_by_sample[sample] = metadata

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
    raw_read_pairs = parse_cutadapt_qc(paths.get("cutadapt_qc", []))
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

    all_chart_samples = list(sample_names)
    for sample in ordered_samples:
        if sample not in all_chart_samples:
            all_chart_samples.append(sample)
    post_dada2_reads = {
        sample: (stats16.get(sample, {}).get("final") or 0)
        + (stats18.get(sample, {}).get("final") or 0)
        for sample in all_chart_samples
        if sample in stats16 or sample in stats18
    }
    raw_reads_chart = svg_sample_read_counts(
        raw_read_pairs, all_chart_samples, "Raw read pairs per sample before filtering and QC"
    )
    post_dada2_chart = svg_sample_read_counts(
        post_dada2_reads, all_chart_samples, "Non-chimeric reads per sample after DADA2"
    )

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
    internal_table_paths = [
        path for path in paths.get("internal_standard_table", []) if Path(path).is_file()
    ]
    interactive_internal = internal_standard_html(
        config, read_tsv(internal_table_paths[0]) if internal_table_paths else []
    )
    internal_section = ""
    if interactive_internal or internal_paths:
        internal_content = interactive_internal or "".join(
            embedded_image(path) for path in internal_paths
        )
        internal_section = f"""
        <section id="internal-standards">
          <div class="eyebrow">Optional analysis</div><h2>Internal-standard correction</h2>
          <p>These offline HTML charts support hover values and horizontal scrolling. Recovery ratios and copy abundances use logarithmic y-axes so large dynamic ranges remain visible. Open a correction method below to inspect its domain-abundance plot.</p>
          {internal_content}
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
.explorer-controls{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:22px 0 12px}} .explorer-control label{{display:block;margin-bottom:5px;color:#344054;font-size:13px;font-weight:750}} .explorer-control select{{width:100%;padding:10px 12px;border:1px solid #cbd5e1;border-radius:8px;background:white;color:var(--ink);font:inherit}} .explorer-control select:disabled{{background:#f1f5f9;color:#94a3b8}} .small-muted{{color:var(--muted);font-size:13px}} .taxonomy-chart{{overflow-x:auto;margin:20px 0 8px;border:1px solid var(--line);border-radius:10px;background:linear-gradient(to top,#f8fafc 1px,transparent 1px);background-size:100% 25%;padding:20px 18px}} .taxonomy-stacked-plot{{display:flex;align-items:flex-start;gap:5px;height:410px;min-width:max-content}} .taxonomy-sample-column{{position:relative;display:flex;flex-direction:column;justify-content:flex-start;width:34px;height:400px}} .taxonomy-stacked-bar{{display:flex;flex:none;flex-direction:column-reverse;width:100%;height:300px;background:#eef2f6;border-bottom:1px solid #94a3b8}} .taxonomy-segment{{width:100%;min-height:1px}} .taxonomy-sample-label{{position:absolute;top:308px;left:16px;width:115px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;transform:rotate(55deg);transform-origin:left top;color:#475467;font-size:11px}} #taxonomy-explorer-legend{{margin:8px 0 24px;min-width:0;padding-left:0}}
.figure-scroll{{overflow-x:auto;margin:18px 0;border:1px solid var(--line);border-radius:10px;background:white}} .report-figure{{display:block;max-width:none;height:auto;margin:0}} code{{white-space:normal;word-break:break-word}} .note{{background:#eff6ff;border-left:4px solid var(--blue);padding:12px 15px;border-radius:6px;color:#344054}}
details{{border:1px solid var(--line);border-radius:10px;margin:12px 0;padding:0 14px 14px}} summary{{cursor:pointer;font-weight:750;padding:14px 0}}
@media(max-width:650px){{.taxonomy-sample-column{{width:30px}}}}
@media print{{nav{{display:none}}body{{background:white}}section{{box-shadow:none;break-inside:avoid}}}}
</style></head><body>
<header><div class="eyebrow" style="color:#bfdbfe">515Y/926R amplicon workflow</div><h1>{esc(study)}</h1><p>Pipeline summary generated {esc(generated)}</p></header>
<nav><a href="#overview">Overview</a><a href="#parameters">Parameters</a><a href="#quality">DADA2</a><a href="#composition">Composition</a><a href="#taxa">Taxa</a>{'<a href="#internal-standards">Internal standards</a>' if interactive_internal or internal_paths else ''}</nav>
<main>
<section id="overview"><div class="eyebrow">Run at a glance</div><h2>Analysis overview</h2>
<div class="cards"><div class="card"><strong>{len(sample_names):,}</strong><span>configured samples</span></div><div class="card"><strong>{fmt_count(split_total)}</strong><span>reads assigned by 16S/18S split</span></div><div class="card"><strong>{fmt_count(final16 + final18)}</strong><span>non-chimeric reads after DADA2</span></div><div class="card"><strong>{len(asvs):,}</strong><span>observed ASVs</span></div><div class="card"><strong>{fmt_percent(median_retention)}</strong><span>median DADA2 retention</span></div></div>
<p class="note">This is a rapid quality-control summary, not a substitute for inspecting unusual samples, QIIME 2 quality visualizations, or the full result tables.</p></section>
<section id="parameters"><div class="eyebrow">Reproducibility</div><h2>Parameters used</h2><h3>Effective DADA2 settings used</h3><p>This table records the values actually applied by the pipeline. For an older config without a <code>dada2</code> block, the workflow defaults are shown.</p><div class="table-wrap"><table><thead><tr><th>Path</th><th>Parameter</th><th>Value</th><th>What it controls</th></tr></thead><tbody>{dada2_parameter_rows}</tbody></table></div><h3>Complete configuration</h3><div class="table-wrap"><table><tbody>{parameter_rows}</tbody></table></div></section>
<section id="quality"><div class="eyebrow">Read processing and DADA2</div><h2>Reads before filtering and quality control</h2><p>These are raw paired-end records reported by Cutadapt before primer removal or other pipeline filtering. One read pair is counted once so it remains comparable with the amplicon counts retained after DADA2.</p>{raw_reads_chart}<h2>Where reads were lost in DADA2</h2><p>Each loss is shown as a read count and the percentage lost from the immediately preceding stage. Filtering covers DADA2 quality filtering and truncation; denoising applies the learned error model; pair merging applies only to paired 16S reads; and the final loss is chimera removal. The 18S reads were concatenated before entering single-end DADA2, so pair merging is not applicable to that path.</p><div class="table-wrap"><table><thead><tr><th>Path</th><th>DADA2 input</th><th>Filtering loss</th><th>Denoising loss</th><th>Pair-merging loss</th><th>Chimera-removal loss</th><th>Final reads</th><th>Total retention</th></tr></thead><tbody>{dada2_summary_rows}</tbody></table></div><h3>DADA2 losses by sample</h3><p>Use this table to identify whether an individual sample loses most reads during filtering, denoising, paired-read merging, or chimera removal. Values below 40% total retention are highlighted for review; these thresholds are guides rather than automatic pass/fail criteria.</p><div class="table-wrap"><table><thead><tr><th>Sample</th><th>Path</th><th>DADA2 input</th><th>Filtering loss</th><th>Denoising loss</th><th>Pair-merging loss</th><th>Chimera-removal loss</th><th>Final reads</th><th>Total retention</th></tr></thead><tbody>{''.join(dada2_sample_rows)}</tbody></table></div><h2>Reads retained after DADA2</h2><p>Each bar is the sample's combined non-chimeric 16S and 18S abundance after DADA2 filtering, denoising, 16S pair merging, and chimera removal.</p>{post_dada2_chart}</section>
<section id="composition"><div class="eyebrow">Basic bar plots</div><h2>Domain composition by sample</h2><p>Bars show relative abundance from <code>{esc(abundance_column)}</code>. Hover over a segment for its value.</p>{domain_chart}<h3>Sequence assignments</h3><p>This breakdown uses the pipeline's <code>Sequence_Type</code> field and taxonomy labels. The broad 16S total includes prokaryotic, chloroplast, and mitochondrial 16S. The figure itself uses mutually exclusive categories, so each sequence count appears in only one bar.</p><div class="cards"><div class="card"><strong>{fmt_count(total_16s)}</strong><span>total 16S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Eukaryotic 18S'])}</strong><span>eukaryotic 18S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Chloroplast 16S'])}</strong><span>chloroplast 16S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Mitochondrial 16S'])}</strong><span>mitochondrial 16S</span></div><div class="card"><strong>{fmt_count(assignment_totals['Unassigned'])}</strong><span>unassigned</span></div></div>{assignment_chart}<h3>Sequence-assignment counts</h3><div class="table-wrap"><table><thead><tr><th>Assignment</th><th>Sequence abundance</th><th>Share of all assignments</th></tr></thead><tbody>{assignment_rows}</tbody></table></div><p class="note"><strong>Total 16S</strong> is a summary row and overlaps its three 16S subcategories; the remaining rows and the figure are mutually exclusive.</p></section>
<section id="taxa"><div class="eyebrow">Taxonomic summary</div><h2>Interactive taxonomy bar plot</h2><p>This QIIME 2-style view shows each sample as a 100% stacked bar. Choose a taxonomy level, then order the bars by SampleID, Condition, Latitude, Longitude, or Depth. When a metadata variable is selected, you can optionally display only one value. Counts use <code>{esc(abundance_column)}</code>; hover over a colored segment for its taxon, relative abundance, and count.</p><div class="explorer-controls"><div class="explorer-control"><label for="taxonomy-rank">Taxonomy level</label><select id="taxonomy-rank"></select></div><div class="explorer-control"><label for="taxonomy-plot-field">Plot samples by</label><select id="taxonomy-plot-field"></select></div><div class="explorer-control"><label for="taxonomy-metadata-value">Filter plotted value</label><select id="taxonomy-metadata-value"></select></div></div><p id="taxonomy-filter-summary" class="small-muted" aria-live="polite"></p><div id="taxonomy-explorer-chart" class="taxonomy-chart" aria-label="Interactive relative taxonomic abundance by sample"></div><div id="taxonomy-explorer-legend" class="chart-legend" aria-label="Taxonomy legend"></div><h3>Top taxa table</h3><div class="table-wrap"><table><thead><tr><th>#</th><th>Taxon</th><th>Total abundance</th><th>Relative abundance</th><th>Samples detected</th></tr></thead><tbody id="taxonomy-explorer-body"></tbody></table></div><noscript><p class="note">Interactive controls require JavaScript. This static summary uses all samples and the first informative SILVA or PR2 rank.</p>{taxa_chart}<div class="table-wrap"><table><thead><tr><th>#</th><th>Taxon</th><th>Total abundance</th><th>Samples detected</th></tr></thead><tbody>{top_taxa_rows}</tbody></table></div></noscript></section>
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
        "cutadapt_qc": list(snakemake_object.input.cutadapt_qc),
        "long_data": str(snakemake_object.input.long_data),
        "internal_standard_figures": list(
            getattr(snakemake_object.input, "internal_standard_figures", []) or []
        ),
        "internal_standard_table": list(
            getattr(snakemake_object.input, "internal_standard_table", []) or []
        ),
    }
    render_report(dict(snakemake_object.config), paths, str(snakemake_object.output.html))


if "snakemake" in globals():
    run_from_snakemake(snakemake)
