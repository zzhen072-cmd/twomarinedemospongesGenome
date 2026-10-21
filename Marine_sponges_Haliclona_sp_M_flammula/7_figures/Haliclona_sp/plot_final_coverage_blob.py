#!/usr/bin/env python3
"""Draw the LHT length-weighted GC/coverage blob plot from an audited TSV.

Required columns: contig_id, length_bp, gc_fraction, ont_coverage, phylum.
No taxonomy is inferred by this script. It only visualizes the supplied labels.
Example:
  python plot_lht_blob.py --table contigs.tsv --out-prefix output/LHT_blob \
      --title 'LHT assembly' --subtitle 'BLAST-derived phylum assignments'

Dependencies: Python 3.9+, numpy, matplotlib. The input TSV is never modified.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from audit_panel_alignment import require_matplotlib_panel_alignment


REQUIRED = ("contig_id", "length_bp", "gc_fraction", "ont_coverage", "phylum")
STATUS_LABELS = {
    "no-hit", "no-hits", "no hit", "no hits", "unresolved", "unclassified",
    "unknown", "na", "n/a", "none", "not-in-supplied-blast",
    "not in supplied blast", "no-taxid", "no taxid", "unresolved-taxonomy",
    "ambiguous", "ambiguous-phylum", "unassigned", "no-phylum",
}
COLORS = (
    "#4C93C3", "#91BD63", "#DB8C4A", "#AD79AD", "#68B5AD",
    "#D96875", "#B4A052", "#7D88BA", "#75ACCD", "#A78B72",
)
RASTER_DPI = 600


def is_status(label):
    lower = label.casefold().strip()
    return lower in STATUS_LABELS or lower.startswith(("unresolved", "not-in-supplied", "no-hit"))


def read_table(path):
    records = []
    seen = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(REQUIRED) - set(reader.fieldnames or []))
        if missing:
            raise ValueError("Missing required TSV columns: " + ", ".join(missing))
        for line, row in enumerate(reader, start=2):
            ident = row["contig_id"].strip()
            if not ident or ident in seen:
                raise ValueError(f"Line {line}: empty or duplicate contig_id {ident!r}")
            seen.add(ident)
            try:
                raw_length = float(row["length_bp"])
                gc = float(row["gc_fraction"])
                coverage = float(row["ont_coverage"])
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Line {line}: invalid length, GC or coverage") from exc
            if not all(math.isfinite(x) for x in (raw_length, gc, coverage)):
                raise ValueError(f"Line {line}: non-finite numeric value")
            if raw_length <= 0 or not raw_length.is_integer():
                raise ValueError(f"Line {line}: length_bp must be a positive integer")
            if not 0 <= gc <= 1 or coverage < 0:
                raise ValueError(f"Line {line}: GC must be 0–1 and coverage must be >= 0")
            label = (row["phylum"] or "").strip()
            if not label:
                raise ValueError(f"Line {line}: phylum is blank; assign an explicit status first")
            records.append({"contig_id": ident, "length_bp": int(raw_length),
                            "gc_fraction": gc, "ont_coverage": coverage, "phylum": label})
    if not records:
        raise ValueError("Input TSV contains no data rows")
    return records


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def stats(records):
    lengths = [r["length_bp"] for r in records]
    plotted = [r for r in records if r["ont_coverage"] > 0]
    return {
        "contig_count": len(records), "sum_length_bp": sum(lengths), "n50_bp": n50(lengths),
        "plotted_contig_count": len(plotted),
        "plotted_sum_length_bp": sum(r["length_bp"] for r in plotted),
        "plotted_n50_bp": n50([r["length_bp"] for r in plotted]),
        "zero_coverage_contig_count": len(records) - len(plotted),
        "zero_coverage_length_bp": sum(r["length_bp"] for r in records if r["ont_coverage"] == 0),
    }


def format_bp(value):
    if value == 0:
        return "0"
    for size, suffix in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(value) >= size:
            scaled = value / size
            return f"{scaled:.1f}{suffix}" if scaled < 100 else f"{scaled:.0f}{suffix}"
    return str(int(value))


def write_tsv(path, fields, records):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def make_groups(records, top_n):
    by_phylum = defaultdict(list)
    for record in records:
        by_phylum[record["phylum"]].append(record)
    ranked = sorted((label for label in by_phylum if not is_status(label)),
                    key=lambda label: (-sum(r["length_bp"] for r in by_phylum[label]), label))
    kept = ranked[:top_n]
    pooled = ranked[top_n:]
    protected = sorted((label for label in by_phylum if is_status(label)),
                       key=lambda label: (-sum(r["length_bp"] for r in by_phylum[label]), label))
    other_label = f"Other ({len(pooled)} phyla)"
    members = {label: [label] for label in kept + protected}
    groups = {label: by_phylum[label] for label in kept + protected}
    order = kept[:]
    if pooled:
        groups[other_label] = [r for label in pooled for r in by_phylum[label]]
        members[other_label] = pooled
        order.append(other_label)
    order += protected
    colors = {label: COLORS[i % len(COLORS)] for i, label in enumerate(kept)}
    if pooled:
        colors[other_label] = "#BAA896"
    # Stable status colors keep subset/reference figures directly comparable.
    for label in protected:
        lower = label.casefold()
        if lower.startswith(("not-in-supplied", "not in supplied")):
            colors[label] = "#CBCFD4"
        elif lower.startswith(("no-hit", "no hit")):
            colors[label] = "#DDE0E4"
        elif lower.startswith("unresolved"):
            colors[label] = "#A2A8B0"
        else:
            colors[label] = "#8A8F97"
    return groups, order, members, colors, pooled


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--title", default="LHT assembly: GC–coverage blob plot")
    parser.add_argument("--subtitle", default="")
    parser.add_argument("--coverage-label", default="ONT coverage (×)",
                        help="Y-axis label identifying the source of the coverage measurement")
    parser.add_argument("--note", default="",
                        help="Replace the final generic footer note (zero-coverage exclusions remain explicit)")
    parser.add_argument("--top-phyla", type=int, default=10)
    parser.add_argument("--gc-bin-width", type=float, default=0.01)
    parser.add_argument("--coverage-bins", type=int, default=40)
    parser.add_argument("--max-bubble-area", type=float, default=740,
                        help="Largest circle marker area in pt²; all others scale strictly by length")
    parser.add_argument("--dpi", type=int, default=RASTER_DPI)
    args = parser.parse_args()
    if not 1 <= args.top_phyla <= len(COLORS):
        parser.error(f"--top-phyla must be between 1 and {len(COLORS)}")
    if not 0 < args.gc_bin_width <= 1:
        parser.error("--gc-bin-width must be > 0 and <= 1")
    if args.coverage_bins < 2 or args.max_bubble_area <= 0 or args.dpi < 72:
        parser.error("Use coverage-bins >= 2, max-bubble-area > 0 and dpi >= 72")
    try:
        records = read_table(args.table)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter

    positive = [r for r in records if r["ont_coverage"] > 0]
    if not positive:
        parser.error("All coverages are zero; a log-coverage plot cannot be drawn")
    groups, order, members, colors, pooled = make_groups(records, args.top_phyla)
    plotted_groups = {label: [r for r in groups[label] if r["ont_coverage"] > 0] for label in order}
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out = lambda suffix: Path(str(args.out_prefix) + suffix)

    summary = [{"plot_group": "Total", "source_phyla": ";".join(sorted({r["phylum"] for r in records})),
                "color_hex": "", **stats(records)}]
    for label in order:
        summary.append({"plot_group": label, "source_phyla": ";".join(members[label]),
                        "color_hex": colors[label], **stats(groups[label])})
    write_tsv(out(".category_summary.tsv"), list(summary[0]), summary)

    # The same coverage-positive set contributes to both marginal histograms.
    gc_values = np.array([r["gc_fraction"] for r in positive])
    gc_width = args.gc_bin_width
    gc_low = max(0, math.floor((float(gc_values.min()) - gc_width) / gc_width) * gc_width)
    gc_high = min(1, math.ceil((float(gc_values.max()) + gc_width) / gc_width) * gc_width)
    gc_bin_count = max(2, int(round((gc_high - gc_low) / gc_width)))
    gc_edges = np.linspace(gc_low, gc_high, gc_bin_count + 1)
    coverage_values = np.array([r["ont_coverage"] for r in positive])
    log_min, log_max = np.log10([coverage_values.min(), coverage_values.max()])
    margin = max(0.08, (log_max - log_min) * 0.045)
    coverage_edges = np.geomspace(10 ** (log_min - margin), 10 ** (log_max + margin), args.coverage_bins + 1)

    hist_gc, hist_coverage, hist_rows = {}, {}, []
    for label in order:
        subset = plotted_groups[label]
        lengths = np.array([r["length_bp"] for r in subset], dtype=np.int64)
        for axis, field, edges, hist_dict in (
            ("gc_fraction", "gc_fraction", gc_edges, hist_gc),
            ("ont_coverage", "ont_coverage", coverage_edges, hist_coverage),
        ):
            values = np.array([r[field] for r in subset])
            weights, _ = np.histogram(values, bins=edges, weights=lengths)
            counts, _ = np.histogram(values, bins=edges)
            hist_dict[label] = weights
            if int(weights.sum()) != int(lengths.sum()):
                raise RuntimeError(f"Histogram lost sequence length: {axis}, {label}")
            for i, (left, right, count, length_sum) in enumerate(zip(edges[:-1], edges[1:], counts, weights)):
                hist_rows.append({"axis": axis, "plot_group": label, "bin_index": i,
                                  "bin_left": f"{left:.15g}", "bin_right": f"{right:.15g}",
                                  "right_edge_inclusive": i == len(edges) - 2,
                                  "contig_count": int(count), "sum_length_bp": int(length_sum)})
    write_tsv(out(".histogram_bins.tsv"), list(hist_rows[0]), hist_rows)

    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 6.0, "axes.linewidth": 0.6,
        "axes.edgecolor": "#6F747B", "axes.labelcolor": "#25282D",
        "xtick.color": "#474D55", "ytick.color": "#474D55", "svg.fonttype": "none",
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
    })
    # Final double-column width: 7.2 in = 182.9 mm.
    fig = plt.figure(figsize=(7.2, 4.65))
    gs = fig.add_gridspec(2, 3, width_ratios=(6.5, 1.65, 4.75), height_ratios=(1.7, 6.5),
                          left=0.073, right=0.984, bottom=0.17, top=0.87,
                          hspace=0.025, wspace=0.045)
    ax = fig.add_subplot(gs[1, 0])
    top = fig.add_subplot(gs[0, 0], sharex=ax)
    right = fig.add_subplot(gs[1, 1], sharey=ax)
    legend_ax = fig.add_subplot(gs[:, 2])
    legend_ax.set_axis_off()
    ax.set_yscale("log")
    ax.set_xlim(gc_edges[0], gc_edges[-1])
    ax.set_ylim(coverage_edges[0], coverage_edges[-1])
    ax.set_xlabel("GC fraction", fontsize=7.2, labelpad=4.5)
    ax.set_ylabel(args.coverage_label, fontsize=7.2, labelpad=5)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.2f}"))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=10))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=tuple(range(2, 10)), numticks=100))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="major", length=4, width=0.8)
    ax.tick_params(which="minor", length=2, width=0.5)
    ax.grid(axis="y", which="major", color="#E7E9ED", linewidth=0.55, zorder=0)
    ax.set_axisbelow(True)
    # One global area multiplier; no minimum-size clamp and no length clipping.
    area_factor = args.max_bubble_area / max(r["length_bp"] for r in positive)
    assignment = {r["contig_id"]: label for label in order for r in plotted_groups[label]}
    draw_records = sorted(positive, key=lambda r: r["length_bp"], reverse=True)
    ax.scatter([r["gc_fraction"] for r in draw_records], [r["ont_coverage"] for r in draw_records],
               s=[r["length_bp"] * area_factor for r in draw_records],
               c=[colors[assignment[r["contig_id"]]] for r in draw_records],
               alpha=0.64, edgecolors="#515B66", linewidths=0.35, zorder=3)

    bottom_gc = np.zeros(len(gc_edges) - 1)
    bottom_cov = np.zeros(len(coverage_edges) - 1)
    for label in order:
        top.bar(gc_edges[:-1], hist_gc[label], width=np.diff(gc_edges), bottom=bottom_gc,
                align="edge", color=colors[label], edgecolor="white", linewidth=0.22, alpha=0.91)
        right.barh(coverage_edges[:-1], hist_coverage[label], height=np.diff(coverage_edges), left=bottom_cov,
                   align="edge", color=colors[label], edgecolor="white", linewidth=0.22, alpha=0.91)
        bottom_gc += hist_gc[label]
        bottom_cov += hist_coverage[label]
    # Both marginal count axes share the same length scale.
    hist_max = max(float(bottom_gc.max()), float(bottom_cov.max())) * 1.10
    top.set_ylim(0, hist_max)
    right.set_xlim(0, hist_max)
    bp_formatter = FuncFormatter(lambda value, _: format_bp(value))
    top.yaxis.set_major_formatter(bp_formatter)
    right.xaxis.set_major_formatter(bp_formatter)
    top.yaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=2))
    right.xaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=2))
    top.set_ylabel("Sum length", fontsize=6.2, labelpad=4)
    right.set_xlabel("Sum length", fontsize=6.2, labelpad=5)
    top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    right.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
    right.tick_params(axis="x", labelsize=5.0, labelrotation=55)
    top.spines[["top", "right"]].set_visible(False)
    right.spines[["top", "right"]].set_visible(False)
    ax.spines[["top", "right"]].set_visible(False)

    # A separate legend panel ensures that long taxonomy labels never overlap data.
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.text(0.05, 0.965, "Phylum / assignment status", fontsize=6.3, weight="bold", va="top")
    col_name, col_count, col_length, col_n50 = 0.092, 0.665, 0.835, 0.992
    for x, label in ((col_count, "Count"), (col_length, "Length"), (col_n50, "N50")):
        legend_ax.text(x, 0.909, label, ha="right", fontsize=5.1, color="#5F6670")
    legend_ax.plot([0.05, 0.995], [0.89, 0.89], color="#CCD0D5", lw=0.7)
    rows = [("Total", positive, "#5F6873")] + [(label, plotted_groups[label], colors[label]) for label in order]
    wrapped = [textwrap.wrap(label, width=24, break_long_words=False, break_on_hyphens=True) for label, _, _ in rows]
    # Allocate line height dynamically, leaving space for the bubble key.
    total_lines = sum(max(1, len(lines)) for lines in wrapped)
    row_unit = min(0.041, 0.62 / max(1, total_lines + len(rows) * 0.25))
    cursor = 0.858
    for (label, subset, color), lines in zip(rows, wrapped):
        lines = lines or [label]
        center_y = cursor - row_unit * (len(lines) - 1) / 2
        legend_ax.scatter([0.06], [center_y], s=14, marker="s", facecolors=color,
                           edgecolors="#616771", linewidths=0.45)
        legend_ax.text(col_name, center_y, "\n".join(lines), ha="left", va="center", fontsize=5.1,
                       linespacing=1.15, weight="bold" if label == "Total" else "normal")
        for x, value in ((col_count, str(len(subset))),
                         (col_length, format_bp(sum(r["length_bp"] for r in subset))),
                         (col_n50, format_bp(n50([r["length_bp"] for r in subset])))):
            legend_ax.text(x, center_y, value, ha="right", va="center", fontsize=5.1,
                           weight="bold" if label == "Total" else "normal")
        cursor -= row_unit * (len(lines) + 0.32)

    key_y = min(0.22, cursor - 0.08)
    legend_ax.text(0.05, key_y + 0.059, "Circle area scales with contig length", fontsize=5.5, weight="bold")
    largest = max(r["length_bp"] for r in positive)
    key_lengths = [largest * ratio for ratio in (0.01, 0.10, 0.50)]
    for x, length in zip((0.18, 0.46, 0.81), key_lengths):
        legend_ax.scatter([x], [key_y - 0.018], s=length * area_factor, color="#BCC2CA",
                           edgecolors="#6F7781", linewidths=0.5, alpha=0.64)
        legend_ax.text(x, key_y - 0.102, format_bp(length) + " bp", ha="center", va="top", fontsize=5.1)

    fig.text(0.073, 0.956, args.title, fontsize=10.5, weight="bold", color="#242B34", va="top")
    if args.subtitle:
        fig.text(0.073, 0.919, "\n".join(textwrap.wrap(args.subtitle, width=150)),
                 fontsize=5.8, color="#59636F", va="top", linespacing=1.25)
    total_stats = stats(records)
    notes = ["Each circle represents one contig; circle area is proportional to sequence length.",
             "Marginal histograms are stacked sums of contig length; coverage uses a logarithmic scale."]
    if pooled:
        notes.append(f"The {args.top_phyla} largest phyla by total sequence length are shown separately; "
                     f"{len(pooled)} additional phyla are combined as Other.")
    zeros = total_stats["zero_coverage_contig_count"]
    if zeros:
        notes.append(f"{zeros:,} zero-coverage contigs ({format_bp(total_stats['zero_coverage_length_bp'])} bp) "
                     "are excluded from this log-scale figure; the summary TSV retains them.")
    if args.note:
        notes.append(args.note)
    elif not zeros:
        notes.append("All table entries are included in the plot and legend statistics.")
    fig.text(0.073, 0.072, "\n".join(notes), fontsize=5.2, color="#59636F", va="top", linespacing=1.35)
    # This is one composite scientific panel: the top/right axes are linked
    # marginals and the final axis is a legend/table, not comparable panels.
    # Audit the single authoritative blob plot area and preserve the N/A report.
    require_matplotlib_panel_alignment(
        fig,
        axes=[ax],
        panel_ids=["blob"],
        json_out=out(".alignment.json"),
        overlay_svg=out(".alignment-overlay.svg"),
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=False,
        strict=True,
    )
    fig.savefig(out(".png"), dpi=args.dpi, facecolor="white")
    fig.savefig(out(".svg"), dpi=args.dpi, facecolor="white")
    fig.savefig(out(".pdf"), dpi=args.dpi, facecolor="white")
    fig.savefig(out(".tiff"), dpi=args.dpi, facecolor="white")
    plt.close(fig)

    # This audit table preserves the exact point sizes and group mapping used.
    point_rows = []
    all_assignment = {r["contig_id"]: label for label in order for r in groups[label]}
    for r in records:
        point_rows.append({**r, "plot_group": all_assignment[r["contig_id"]],
                           "color_hex": colors[all_assignment[r["contig_id"]]],
                           "plotted": r["ont_coverage"] > 0,
                           "marker_area_pt2": f"{r['length_bp'] * area_factor:.12g}"})
    write_tsv(out(".plotted_points.tsv"), list(point_rows[0]), point_rows)
    print(f"Rendered {len(positive):,} of {len(records):,} contigs; {len(order)} displayed groups.")
    print(f"Outputs: {args.out_prefix}.{{png,svg,pdf,tiff,category_summary.tsv,histogram_bins.tsv,plotted_points.tsv}}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
