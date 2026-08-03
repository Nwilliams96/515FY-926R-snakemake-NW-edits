"""Generate a self-contained HTML summary for a completed amplicon run."""

from __future__ import annotations

import base64
import csv
import html
import math
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


PALETTE = [
    "#2563eb", "#0d9488", "#f59e0b", "#dc2626", "#7c3aed",
    "#0891b2", "#65a30d", "#db2777", "#4f46e5", "#ea580c",
    "#64748b",
]


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
        final_reads = number(first_value(row, ["non-chimeric", "non_chimeric", "nonchimeric"]))
        if sample and input_reads is not None and final_reads is not None:
            result[normalize_sample_id(sample)] = {
                "input": input_reads,
                "final": final_reads,
                "retention": final_reads / input_reads if input_reads else None,
            }
    return result


def parse_cutadapt_logs(paths):
    total = 0
    retained = 0
    found_total = False
    found_retained = False
    for path in paths:
        try:
            text = Path(path).read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        total_match = re.search(r"Total read pairs processed:\s*([\d,]+)", text, re.I)
        retained_match = re.search(r"Pairs written \(passing filters\):\s*([\d,]+)", text, re.I)
        if total_match:
            total += int(total_match.group(1).replace(",", ""))
            found_total = True
        if retained_match:
            retained += int(retained_match.group(1).replace(",", ""))
            found_retained = True
    return total if found_total else None, retained if found_retained else None


def parse_quality_profiles(directories):
    """Extract median per-base Phred profiles from QIIME 2 demux QZVs."""
    profiles = []
    directory_labels = ["16S paired", "18S paired", "18S concatenated"]
    for directory_index, directory in enumerate(directories):
        label_prefix = directory_labels[directory_index] if directory_index < len(directory_labels) else Path(directory).name
        for qzv_path in sorted(Path(directory).glob("*.qzv")):
            try:
                archive = zipfile.ZipFile(qzv_path)
            except (OSError, zipfile.BadZipFile):
                continue
            with archive:
                summary_names = [
                    name for name in archive.namelist()
                    if name.lower().endswith("seven-number-summaries.tsv")
                ]
                for summary_name in summary_names:
                    direction_name = Path(summary_name).name.split("-seven-number", 1)[0]
                    direction = direction_name.replace("_", " ").replace("-", " ").strip().title()
                    text = archive.read(summary_name).decode("utf-8-sig", errors="replace")
                    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
                    points = []
                    for row in rows:
                        position = number(first_value(row, ["position", "pos", "base position"]))
                        median = number(first_value(row, ["50%", "median", "50th percentile"]))
                        if position is not None and median is not None:
                            points.append((position, median))
                    if points:
                        profiles.append((f"{label_prefix} {direction}".strip(), points))
    return profiles


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


def svg_read_fate(stages):
    width = 940
    height = 280
    left = 185
    plot_width = 700
    bar_height = 38
    gap = 18
    max_value = max((sum(parts.values()) for _, parts in stages), default=1) or 1
    colors = {"16S": "#2563eb", "18S": "#0d9488", "Reads": "#64748b"}
    pieces = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Read retention through the pipeline">']
    for index, (label, parts) in enumerate(stages):
        y = 24 + index * (bar_height + gap)
        total = sum(parts.values())
        pieces.append(f'<text x="{left - 12}" y="{y + 24}" text-anchor="end" class="svg-label">{esc(label)}</text>')
        x = left
        for category, value in parts.items():
            segment = plot_width * value / max_value
            pieces.append(
                f'<rect x="{x:.1f}" y="{y}" width="{segment:.1f}" height="{bar_height}" '
                f'rx="5" fill="{colors.get(category, "#64748b")}"><title>{esc(category)}: {fmt_count(value)}</title></rect>'
            )
            if segment > 74:
                pieces.append(
                    f'<text x="{x + segment / 2:.1f}" y="{y + 24}" text-anchor="middle" class="svg-inside">'
                    f'{esc(category)} {fmt_count(value)}</text>'
                )
            x += segment
        pieces.append(f'<text x="{left + plot_width + 12}" y="{y + 24}" class="svg-total">{fmt_count(total)}</text>')
    pieces.append('</svg>')
    return "".join(pieces)


def svg_quality_profiles(profiles, trunc_r1=None, trunc_r2=None):
    if not profiles:
        return '<p class="note">Median base-quality profiles were not found in the QIIME 2 visualization archives.</p>'
    width, height = 940, 390
    left, top, plot_w, plot_h = 65, 28, 830, 285
    max_x = max(point[0] for _, points in profiles for point in points) or 1
    max_y = max(45, math.ceil(max(point[1] for _, points in profiles for point in points) / 5) * 5)
    pieces = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Median per-base Phred quality">']
    for tick in range(0, int(max_y) + 1, 5):
        y = top + plot_h * (1 - tick / max_y)
        pieces.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>')
        pieces.append(f'<text x="{left - 9}" y="{y + 4:.1f}" text-anchor="end" class="svg-label">{tick}</text>')
    for truncation, label in ((number(trunc_r1), "R1 truncation"), (number(trunc_r2), "R2 truncation")):
        if truncation is not None and truncation <= max_x:
            x = left + plot_w * truncation / max_x
            pieces.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" class="trunc-line"/>')
            pieces.append(f'<text x="{x + 4:.1f}" y="{top + 13}" class="trunc-label">{esc(label)}</text>')
    for index, (label, points) in enumerate(profiles):
        color = PALETTE[index % len(PALETTE)]
        coordinates = " ".join(
            f"{left + plot_w * x / max_x:.1f},{top + plot_h * (1 - y / max_y):.1f}"
            for x, y in points
        )
        pieces.append(f'<polyline points="{coordinates}" fill="none" stroke="{color}" stroke-width="2.2"><title>{esc(label)}</title></polyline>')
        legend_x = left + (index % 3) * 270
        legend_y = 345 + (index // 3) * 22
        pieces.append(f'<line x1="{legend_x}" y1="{legend_y - 4}" x2="{legend_x + 24}" y2="{legend_y - 4}" stroke="{color}" stroke-width="3"/>')
        pieces.append(f'<text x="{legend_x + 31}" y="{legend_y}" class="svg-label">{esc(label)}</text>')
    pieces.append(f'<text x="{left + plot_w / 2}" y="{top + plot_h + 31}" text-anchor="middle" class="svg-label">Base position</text>')
    pieces.append(f'<text transform="translate(17,{top + plot_h / 2}) rotate(-90)" text-anchor="middle" class="svg-label">Median Phred score</text>')
    pieces.append('</svg>')
    return "".join(pieces)


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
    height = 430
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
    legend_y = 390
    legend_x = 65
    for index, category in enumerate(categories):
        x = legend_x + (index % 5) * 145
        y = legend_y + (index // 5) * 24
        pieces.append(f'<rect x="{x}" y="{y - 11}" width="12" height="12" fill="{PALETTE[index % len(PALETTE)]}"/>')
        pieces.append(f'<text x="{x + 18}" y="{y}" class="svg-label">{esc(category)}</text>')
    pieces.append('</svg></div>')
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


def quality_status(retention):
    if retention is None:
        return "not available", "neutral"
    if retention >= 0.70:
        return "good", "good"
    if retention >= 0.40:
        return "review", "warn"
    return "low", "bad"


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
    raw_pairs, primer_retained = parse_cutadapt_logs(paths.get("trimming_logs", []))
    quality_profiles = parse_quality_profiles(paths.get("quality_directories", []))

    long_rows = read_tsv(paths["long_data"])
    abundance_column = "Corrected_Sequence_Counts"
    if long_rows and abundance_column not in long_rows[0]:
        abundance_column = "Raw_Sequence_Counts"
    domains_by_sample = defaultdict(Counter)
    taxa_totals = Counter()
    taxa_samples = defaultdict(set)
    asvs = set()
    for row in long_rows:
        sample = normalize_sample_id(row.get("SampleID"))
        abundance = number(row.get(abundance_column)) or 0
        if not sample or abundance <= 0:
            continue
        domain = (row.get("Domain") or "Unassigned").strip() or "Unassigned"
        taxon = taxonomy_label(row)
        domains_by_sample[sample][domain] += abundance
        taxa_totals[taxon] += abundance
        taxa_samples[taxon].add(sample)
        if row.get("ASV_hash"):
            asvs.add(row["ASV_hash"])

    ordered_samples = [sample for sample in sample_names if sample in set(split) | set(stats16) | set(stats18)]
    for sample in sorted((set(split) | set(stats16) | set(stats18)) - set(ordered_samples)):
        ordered_samples.append(sample)

    split_total = sum(item["total"] for item in split.values())
    prok_split = sum(item["prok"] for item in split.values())
    euk_split = sum(item["euk"] for item in split.values())
    final16 = sum(item["final"] for item in stats16.values())
    final18 = sum(item["final"] for item in stats18.values())
    all_retentions = [item["retention"] for item in list(stats16.values()) + list(stats18.values()) if item["retention"] is not None]
    median_retention = None
    if all_retentions:
        ordered = sorted(all_retentions)
        middle = len(ordered) // 2
        median_retention = ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2

    stages = []
    if raw_pairs is not None:
        stages.append(("Raw paired reads", {"Reads": raw_pairs}))
    if primer_retained is not None:
        stages.append(("Primer-matched reads", {"Reads": primer_retained}))
    stages.extend([
        ("16S / 18S split", {"16S": prok_split, "18S": euk_split}),
        ("After DADA2", {"16S": final16, "18S": final18}),
    ])

    quality_rows = []
    for sample in ordered_samples:
        split_row = split.get(sample, {})
        s16 = stats16.get(sample, {})
        s18 = stats18.get(sample, {})
        status16, class16 = quality_status(s16.get("retention"))
        status18, class18 = quality_status(s18.get("retention"))
        total = split_row.get("total")
        euk_fraction = split_row.get("euk", 0) / total if total else None
        quality_rows.append(
            f"<tr><td>{esc(sample)}</td><td>{fmt_count(total)}</td><td>{fmt_percent(euk_fraction)}</td>"
            f"<td>{fmt_count(s16.get('input'))}</td><td>{fmt_count(s16.get('final'))}</td>"
            f"<td><span class='pill {class16}'>{fmt_percent(s16.get('retention'))} · {status16}</span></td>"
            f"<td>{fmt_count(s18.get('input'))}</td><td>{fmt_count(s18.get('final'))}</td>"
            f"<td><span class='pill {class18}'>{fmt_percent(s18.get('retention'))} · {status18}</span></td></tr>"
        )

    parameter_rows = "".join(
        f"<tr><th>{esc(key)}</th><td><code>{esc(value)}</code></td></tr>"
        for key, value in flatten_config(config)
    )
    top_taxa_rows = "".join(
        f"<tr><td>{index}</td><td>{esc(taxon)}</td><td>{fmt_count(value)}</td><td>{len(taxa_samples[taxon])}</td></tr>"
        for index, (taxon, value) in enumerate(taxa_totals.most_common(20), 1)
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
    taxa_chart = svg_top_taxa(taxa_totals)
    truncation_config = config.get("trunclens", {})
    quality_chart = svg_quality_profiles(
        quality_profiles,
        truncation_config.get("truncR1"),
        truncation_config.get("truncR2"),
    )
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
.svg-label{{fill:#475467;font-size:12px}} .svg-total{{fill:#344054;font-size:12px;font-weight:700}} .svg-inside{{fill:white;font-size:11px;font-weight:700}} .grid{{stroke:#e4e7ec;stroke-width:1}} svg{{max-width:100%;height:auto}} .chart-scroll{{overflow-x:auto}}
.trunc-line{{stroke:#b42318;stroke-width:1.5;stroke-dasharray:5 4}} .trunc-label{{fill:#b42318;font-size:10px;font-weight:700}}
.report-figure{{display:block;max-width:100%;height:auto;border:1px solid var(--line);border-radius:10px;margin:18px auto}} code{{white-space:normal;word-break:break-word}} .note{{background:#eff6ff;border-left:4px solid var(--blue);padding:12px 15px;border-radius:6px;color:#344054}}
@media print{{nav{{display:none}}body{{background:white}}section{{box-shadow:none;break-inside:avoid}}}}
</style></head><body>
<header><div class="eyebrow" style="color:#bfdbfe">515Y/926R amplicon workflow</div><h1>{esc(study)}</h1><p>Pipeline summary generated {esc(generated)}</p></header>
<nav><a href="#overview">Overview</a><a href="#parameters">Parameters</a><a href="#read-fate">Read fate</a><a href="#quality">Quality</a><a href="#composition">Composition</a><a href="#taxa">Taxa</a>{'<a href="#internal-standards">Internal standards</a>' if internal_paths else ''}</nav>
<main>
<section id="overview"><div class="eyebrow">Run at a glance</div><h2>Analysis overview</h2>
<div class="cards"><div class="card"><strong>{len(sample_names):,}</strong><span>configured samples</span></div><div class="card"><strong>{fmt_count(split_total)}</strong><span>reads assigned by 16S/18S split</span></div><div class="card"><strong>{fmt_count(final16 + final18)}</strong><span>non-chimeric reads after DADA2</span></div><div class="card"><strong>{len(asvs):,}</strong><span>observed ASVs</span></div><div class="card"><strong>{fmt_percent(median_retention)}</strong><span>median DADA2 retention</span></div></div>
<p class="note">This is a rapid quality-control summary, not a substitute for inspecting unusual samples, QIIME 2 quality visualizations, or the full result tables.</p></section>
<section id="parameters"><div class="eyebrow">Reproducibility</div><h2>Parameters used</h2><div class="table-wrap"><table><tbody>{parameter_rows}</tbody></table></div></section>
<section id="read-fate"><div class="eyebrow">Processing losses</div><h2>Where reads were retained or lost</h2><p>Counts are aggregated across samples. The 16S and 18S paths branch after database-based read splitting.</p>{svg_read_fate(stages)}</section>
<section id="quality"><div class="eyebrow">Per-base and per-sample QC</div><h2>Read quality and DADA2 retention</h2><h3>Median base-quality profiles</h3><p>Profiles are extracted from the QIIME 2 demultiplexing visualizations. Dashed lines show the configured DADA2 truncation positions.</p>{quality_chart}<h3>DADA2 retention by sample</h3><p>Retention is the non-chimeric read count divided by the DADA2 input count. Values below 40% are highlighted for review; thresholds are guides rather than pass/fail criteria.</p><div class="table-wrap"><table><thead><tr><th>Sample</th><th>Split reads</th><th>18S fraction</th><th>16S input</th><th>16S final</th><th>16S retention</th><th>18S input</th><th>18S final</th><th>18S retention</th></tr></thead><tbody>{''.join(quality_rows)}</tbody></table></div></section>
<section id="composition"><div class="eyebrow">Basic bar plot</div><h2>Domain composition by sample</h2><p>Bars show relative abundance from <code>{esc(abundance_column)}</code>. Hover over a segment for its value.</p>{domain_chart}</section>
<section id="taxa"><div class="eyebrow">Taxonomic summary</div><h2>Most abundant taxa</h2><p>For SILVA assignments this uses phylum where available; for PR2 it uses division or supergroup.</p>{taxa_chart}<h3>Top taxa table</h3><div class="table-wrap"><table><thead><tr><th>Rank</th><th>Taxon</th><th>Total abundance</th><th>Samples detected</th></tr></thead><tbody>{top_taxa_rows}</tbody></table></div></section>
{internal_section}
</main></body></html>"""
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
        "quality_directories": list(
            getattr(snakemake_object.input, "quality_directories", []) or []
        ),
        "trimming_logs": list(getattr(snakemake_object.input, "trimming_logs", []) or []),
        "internal_standard_figures": list(
            getattr(snakemake_object.input, "internal_standard_figures", []) or []
        ),
    }
    render_report(dict(snakemake_object.config), paths, str(snakemake_object.output.html))


if "snakemake" in globals():
    run_from_snakemake(snakemake)
