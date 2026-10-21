#!/usr/bin/env python3
"""Plot verified gene-level annotation summaries for LHT-04 and NJ-lan.

No illustrative data are generated. All four required TSVs must exist. UpSet
sources are defined by ``source_order`` in upset_intersections.tsv, not inferred
from the coverage table (which may contain additional annotation measures).

Example:
  python plot_functional_main.py --tables ../tables --out ../figures
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
import textwrap
from pathlib import Path

from audit_panel_alignment import require_matplotlib_panel_alignment


SPECIES = ("LHT-04", "NJ-lan")
COLORS = {"LHT-04": "#277D91", "NJ-lan": "#D18136"}
LABELS = {"LHT-04": "LHT-04", "NJ-lan": "NJ-lan"}
ONTOLOGIES = ("BP", "MF", "CC")
ONTOLOGY_NAMES = {
    "BP": "Biological process", "MF": "Molecular function",
    "CC": "Cellular component",
}
DISPLAY_SOURCES = {"SwissProt": "Swiss-Prot"}
REQUIRED = {
    "functional_coverage.tsv": {"species", "source", "gene_count", "total_genes", "percent_total_genes"},
    "upset_intersections.tsv": {"species", "pattern", "sources", "source_order", "gene_count", "total_genes", "percent_total_genes"},
    "goslim_counts.tsv": {"species", "ontology", "GO_ID", "GO_term", "gene_count", "total_genes", "percent_total_genes"},
    "kegg_level2_counts.tsv": {"species", "A_id", "A_name", "B_id", "B_name", "gene_count", "total_genes", "percent_total_genes"},
}


class DataError(ValueError):
    """A missing or internally inconsistent input table."""


def read_table(path: Path, required: set[str]) -> list[dict]:
    if not path.is_file():
        raise DataError(f"Missing required input: {path}. No figures were drawn.")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        absent = required - set(reader.fieldnames or [])
        if absent:
            raise DataError(f"{path.name}: missing columns: {', '.join(sorted(absent))}")
        rows = list(reader)
    if not rows:
        raise DataError(f"{path.name}: contains no data rows; cannot draw a verified figure.")
    for line, row in enumerate(rows, 2):
        if row["species"] not in SPECIES:
            raise DataError(f"{path.name}:{line}: unexpected species {row['species']!r}")
        try:
            row["gene_count"] = int(row["gene_count"])
            row["total_genes"] = int(row["total_genes"])
            row["percent_total_genes"] = float(row["percent_total_genes"])
        except (ValueError, TypeError) as exc:
            raise DataError(f"{path.name}:{line}: non-numeric count or percentage") from exc
        n, total, pct = (row[k] for k in ("gene_count", "total_genes", "percent_total_genes"))
        if total <= 0 or not 0 <= n <= total or not math.isfinite(pct):
            raise DataError(f"{path.name}:{line}: invalid count, denominator or percentage")
        if abs(pct - 100 * n / total) > 0.011:
            raise DataError(f"{path.name}:{line}: percentage differs from 100 × count / total")
    return rows


def keyed(rows: list[dict], fields: tuple[str, ...], name: str) -> dict:
    output = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        if key in output:
            raise DataError(f"{name}: duplicated row key {key}")
        output[key] = row
    return output


def normalize_ontology(value: str) -> str:
    value = value.strip().lower().replace("_", " ")
    aliases = {"bp": "BP", "biological process": "BP", "mf": "MF", "molecular function": "MF",
               "cc": "CC", "cellular component": "CC"}
    if value not in aliases:
        raise DataError(f"Unknown GO ontology {value!r}; expected BP/MF/CC")
    return aliases[value]


def load_data(tables: Path) -> dict:
    # Validate every required file before importing the plotting library or
    # creating output directories, so absent data cannot yield partial figures.
    missing = [str(tables / name) for name in REQUIRED if not (tables / name).is_file()]
    if missing:
        raise DataError("Missing required annotation tables:\n  " + "\n  ".join(missing)
                        + "\nNo figures were drawn; real gene-level annotation input is required.")
    data = {name: read_table(tables / name, fields) for name, fields in REQUIRED.items()}
    totals = {}
    for name, rows in data.items():
        if {r["species"] for r in rows} != set(SPECIES):
            raise DataError(f"{name}: must contain both LHT-04 and NJ-lan")
        for row in rows:
            sp, total = row["species"], row["total_genes"]
            if sp in totals and totals[sp] != total:
                raise DataError(f"Inconsistent total_genes for {sp}: {totals[sp]} versus {total} in {name}")
            totals[sp] = total
    overlap = data["upset_intersections.tsv"]
    orders = {r["source_order"] for r in overlap}
    if len(orders) != 1:
        raise DataError("upset_intersections.tsv: source_order must be identical for both species and every row")
    sources = tuple(next(iter(orders)).split(";"))
    if sources not in (("eggNOG", "InterPro"), ("eggNOG", "InterPro", "SwissProt")):
        raise DataError("UpSet source_order must be eggNOG;InterPro or eggNOG;InterPro;SwissProt")
    lookup = keyed(overlap, ("species", "pattern"), "upset_intersections.tsv")
    expected_patterns = {"".join(bits) for bits in itertools.product("01", repeat=len(sources))}
    for sp in SPECIES:
        observed = {r["pattern"] for r in overlap if r["species"] == sp}
        if observed != expected_patterns:
            raise DataError(f"UpSet {sp}: expected all {2 ** len(sources)} binary patterns, including none; "
                            f"missing={sorted(expected_patterns - observed)}, extra={sorted(observed - expected_patterns)}")
        if sum(lookup[(sp, p)]["gene_count"] for p in expected_patterns) != totals[sp]:
            raise DataError(f"UpSet {sp}: exclusive intersections do not sum to total_genes")
    for row in overlap:
        active = tuple(s for bit, s in zip(row["pattern"], sources) if bit == "1")
        supplied = tuple(s.strip() for s in row["sources"].split(";") if s.strip())
        if not active:
            if supplied not in ((), ("None",), ("none",)):
                raise DataError(f"UpSet {row['species']}: none pattern has sources={row['sources']!r}")
        elif supplied != active:
            raise DataError(f"UpSet sources do not match pattern={row['pattern']}: {row['sources']!r}")
    coverage = keyed(data["functional_coverage.tsv"], ("species", "source"), "functional_coverage.tsv")
    set_sizes = {}
    for sp in SPECIES:
        for i, source in enumerate(sources):
            marginal = sum(r["gene_count"] for r in overlap if r["species"] == sp and r["pattern"][i] == "1")
            set_sizes[(sp, source)] = marginal
            if (sp, source) in coverage and coverage[(sp, source)]["gene_count"] != marginal:
                raise DataError(f"Coverage/UpSet mismatch for {sp} {source}")
    set_path = tables / "upset_set_sizes.tsv"
    if set_path.is_file():
        for row in read_table(set_path, REQUIRED["functional_coverage.tsv"]):
            key = (row["species"], row["source"])
            if key not in set_sizes or set_sizes[key] != row["gene_count"] or totals[row["species"]] != row["total_genes"]:
                raise DataError(f"upset_set_sizes.tsv does not agree with intersections: {key}")
    go = data["goslim_counts.tsv"]
    for row in go:
        row["ontology"] = normalize_ontology(row["ontology"])
    keyed(go, ("species", "GO_ID"), "goslim_counts.tsv")
    keyed(data["kegg_level2_counts.tsv"], ("species", "A_id", "B_id"), "kegg_level2_counts.tsv")
    # A category may be absent from one species (treated as zero for plotting),
    # but identity, label and ontology must agree when present in both species.
    for rows, keys, labels in ((go, ("GO_ID",), ("ontology", "GO_term")),
                               (data["kegg_level2_counts.tsv"], ("A_id", "B_id"), ("A_name", "B_name"))):
        names = {}
        for row in rows:
            key, label = tuple(row[k] for k in keys), tuple(row[k] for k in labels)
            if key in names and names[key] != label:
                raise DataError(f"Category identity/label mismatch for {key}")
            names[key] = label
    none = "0" * len(sources)
    patterns = sorted(expected_patterns - {none},
                      key=lambda p: (-sum(lookup[(sp, p)]["percent_total_genes"] for sp in SPECIES), -p.count("1"), p))
    patterns.append(none)
    return {"tables": data, "totals": totals, "sources": sources, "overlap": lookup,
            "patterns": patterns, "set_sizes": set_sizes}


def category_records(rows: list[dict], keys: tuple[str, ...], label: str, group: str) -> list[dict]:
    records = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        record = records.setdefault(key, {"key": key, "label": row[label], "group": row[group],
                                          "values": {sp: 0.0 for sp in SPECIES}, "counts": {sp: 0 for sp in SPECIES}})
        record["values"][row["species"]] = row["percent_total_genes"]
        record["counts"][row["species"]] = row["gene_count"]
    for record in records.values():
        record["mean"] = sum(record["values"].values()) / len(SPECIES)
    return list(records.values())


def configure_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib import font_manager
    installed = {font.name for font in font_manager.fontManager.ttflist}
    family = "Arial" if "Arial" in installed else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family, "font.size": 10, "axes.labelsize": 9,
        "axes.titlesize": 11, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    return plt, family


def panel_letter(fig, rect, letter):
    fig.text(rect[0] - 0.035, rect[1] + rect[3] + 0.023, letter,
             fontsize=14, fontweight="bold", ha="left", va="bottom")


def draw_upset(fig, rect, sp, data, ymax, letter):
    """Compact UpSet with marginal set sizes in row labels and a none column."""
    left, bottom, width, height = rect
    patterns, sources = data["patterns"], data["sources"]
    matrix_height = 0.073 if len(sources) == 2 else 0.090
    matrix = fig.add_axes([left, bottom, width, matrix_height])
    bars = fig.add_axes([left, bottom + matrix_height + 0.031, width,
                         height - matrix_height - 0.031])
    color, total = COLORS[sp], data["totals"][sp]
    values = [data["overlap"][(sp, p)]["percent_total_genes"] for p in patterns]
    counts = [data["overlap"][(sp, p)]["gene_count"] for p in patterns]
    bars.bar(range(len(patterns)), values,
             color=[color if "1" in p else "#B6BCC1" for p in patterns], width=0.68, zorder=3)
    bars.set_xlim(-.6, len(patterns) - .4)
    bars.set_ylim(0, ymax)
    bars.set_ylabel("All predicted genes (%)", labelpad=5)
    bars.set_xticks([])
    bars.grid(False)
    bars.spines["bottom"].set_visible(False)
    bars.set_title(LABELS[sp], loc="left", fontweight="bold", color=color, pad=26)
    bars.text(0, 1.035, f"{total:,} representative proteins",
              transform=bars.transAxes, fontsize=8.5, color="#515B63", va="bottom")
    for j in range(len(sources)):
        if j % 2 == 0:
            matrix.axhspan(j - .43, j + .43, color="#F3F5F6", zorder=0)
    for i, pattern in enumerate(patterns):
        active = [j for j, bit in enumerate(pattern) if bit == "1"]
        matrix.scatter([i] * len(sources), range(len(sources)), s=21, color="#D9DEE1", zorder=2)
        if active:
            matrix.plot([i, i], [min(active), max(active)], color=color, lw=1.6, zorder=3)
            matrix.scatter([i] * len(active), active, s=27, color=color, zorder=4)
    matrix.set_xlim(bars.get_xlim())
    matrix.set_ylim(len(sources) - .45, -.55)
    matrix.set_yticks(range(len(sources)))
    matrix.set_yticklabels([
        f"{DISPLAY_SOURCES.get(s, s)}  {100 * data['set_sizes'][(sp, s)] / total:.1f}%"
        for s in sources
    ], fontsize=8)
    matrix.set_xticks(range(len(patterns)))
    matrix.set_xticklabels(["None" if "1" not in p else "" for p in patterns], fontsize=8)
    matrix.tick_params(axis="both", length=0, pad=3)
    for spine in matrix.spines.values():
        spine.set_visible(False)
    bars.text(-.32, 1.19, letter, transform=bars.transAxes, fontsize=14,
              fontweight="bold", ha="left", va="bottom", clip_on=False)
    return bars


def draw_categories(ax, records, title, group_names=None, wrap=27, group_headers=True):
    """Paired points using the same category set and denominator for each sample."""
    positions, labels, headers = [], [], []
    cursor = 0.0
    prior_group = None
    for record in records:
        if group_headers and record["group"] != prior_group:
            if prior_group is not None:
                cursor += .55
            headers.append((cursor, (group_names or {}).get(record["group"], record["group"])))
            cursor += .70
            prior_group = record["group"]
        wrapped = textwrap.fill(record["label"], wrap, break_long_words=False, break_on_hyphens=False)
        positions.append(cursor)
        labels.append(wrapped)
        cursor += max(1.05, .54 * (wrapped.count("\n") + 1))
    for i, (pos, record) in enumerate(zip(positions, records)):
        if i % 2 == 0:
            ax.axhspan(pos - .45, pos + .45, color="#F5F7F8", zorder=0)
        first, second = (record["values"][sp] for sp in SPECIES)
        ax.plot([first, second], [pos - .11, pos + .11], color="#AEB8BE", lw=1.05, zorder=2)
        for sp, offset in (("LHT-04", -.11), ("NJ-lan", .11)):
            ax.scatter(record["values"][sp], pos + offset, s=25, color=COLORS[sp],
                       edgecolors="white", linewidths=.45, zorder=3)
    max_value = max((max(r["values"].values()) for r in records), default=0)
    ax.set_xlim(-max(.4, max_value * .025), max(2, max_value * 1.13))
    ax.set_ylim(cursor - .20, -.55)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8.2, linespacing=1.15)
    ax.tick_params(axis="y", length=0, pad=7)
    ax.set_xlabel("All predicted genes (%)", labelpad=7)
    ax.grid(axis="x", color="#E3E7EA", lw=.6, zorder=0)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#AAB1B5")
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", pad=17)
    for pos, label in headers:
        ax.text(-.025, pos - .05, label, transform=ax.get_yaxis_transform(), fontsize=8.4,
                fontweight="bold", color="#46545D", ha="right", va="center", clip_on=False)
    return positions


def add_legend(fig, y):
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", markersize=5.5, linestyle="none", color=COLORS[sp], label=LABELS[sp])
               for sp in SPECIES]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.5, y), frameon=False,
               ncol=2, fontsize=9, handletextpad=.4, columnspacing=2.8)


def save_figure(fig, path: Path, plt, dpi: int) -> list[str]:
    minimum_dpi = 600
    dpi = max(dpi, minimum_dpi)
    outputs = [path.with_suffix(".png"), path.with_suffix(".svg"), path.with_suffix(".pdf"), path.with_suffix(".tiff")]
    # Preserve the stated physical figure dimensions; the layout uses
    # explicit margins to keep labels within the canvas.
    fig.savefig(outputs[0], dpi=dpi, facecolor="white")
    fig.savefig(outputs[1], dpi=dpi, facecolor="white")
    fig.savefig(outputs[2], dpi=dpi, facecolor="white")
    fig.savefig(outputs[3], dpi=dpi, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    files = [str(output.resolve()) for output in outputs]
    plt.close(fig)
    return files


def main_figure(data, go, kegg, out, plt, dpi):
    fig = plt.figure(figsize=(7.2, 10.8))
    fig.text(.06, .980, "Functional annotation of the predicted gene sets", fontsize=14,
             fontweight="bold", ha="left", va="top")
    ymax_raw = max(r["percent_total_genes"] for r in data["overlap"].values())
    ymax = math.ceil(max(10, ymax_raw * 1.17) / 10) * 10
    top_rects = {"LHT-04": (.140, .704, .338, .195), "NJ-lan": (.641, .704, .338, .195)}
    upset_axes = []
    for sp, letter in zip(SPECIES, ("a", "b")):
        upset_axes.append(draw_upset(fig, top_rects[sp], sp, data, ymax, letter))
    fig.text(.06, .670, "Source rows: marginal annotation coverage. Grey column: no match in the displayed sources.",
             fontsize=8.1, color="#56616A", ha="left")
    add_legend(fig, .642)
    go_rect = (.267, .148, .215, .437)
    kegg_rect = (.790, .148, .188, .437)
    ax_go, ax_kegg = fig.add_axes(go_rect), fig.add_axes(kegg_rect)
    draw_categories(ax_go, go, "GO slim", ONTOLOGY_NAMES, wrap=27, group_headers=True)
    draw_categories(ax_kegg, kegg, "KEGG level 2", wrap=27, group_headers=False)
    fig.text(.06, .604, "c", fontsize=14, fontweight="bold")
    fig.text(.56, .604, "d", fontsize=14, fontweight="bold")
    fig.text(.06, .080, "One longest protein per gene. GO slim and KEGG categories are non-exclusive; percentages use all genes.",
             fontsize=8.1, color="#4D5961")
    fig.text(.06, .062, "Categories were selected by mean coverage across the two species. Human-disease categories are omitted from panel D.",
             fontsize=8.1, color="#4D5961")
    fig.text(.06, .044, "Descriptive annotation summaries; no enrichment tests. InterProScan denotes any member-database signature.",
             fontsize=8.1, color="#4D5961")
    require_matplotlib_panel_alignment(
        fig,
        axes=upset_axes,
        panel_ids=("a", "b"),
        row_groups=(("a", "b"),),
        json_out=out / "Figure3_functional_annotation.alignment.json",
        overlay_svg=out / "Figure3_functional_annotation.alignment.svg",
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )
    return save_figure(fig, out / "Figure3_functional_annotation", plt, dpi)


def supplement(records, title, path, plt, dpi, group_names=None):
    from matplotlib.backends.backend_pdf import PdfPages
    pages=[];current=[];units=0;last_group=None
    for row in records:
        extra=max(1.05,.54*(textwrap.fill(row['label'],46).count('\n')+1))+(1.25 if row['group']!=last_group else 0)
        if current and units+extra>32:
            pages.append(current);current=[];units=0;last_group=None
            extra=max(1.05,.54*(textwrap.fill(row['label'],46).count('\n')+1))+1.25
        current.append(row);units+=extra;last_group=row['group']
    if current:pages.append(current)
    files=[];maximum=max(max(r['values'].values()) for r in records)
    with PdfPages(path.with_suffix('.pdf')) as pdf:
        for i,page in enumerate(pages,1):
            fig=plt.figure(figsize=(9.4,11.7))
            ax=fig.add_axes([.43,.13,.545,.72])
            draw_categories(ax,page,'',group_names,wrap=46,group_headers=True)
            ax.set_xlim(-max(.4,maximum*.025),max(2,maximum*1.13))
            fig.text(.045,.965,title,fontsize=13,weight='bold',va='top')
            fig.text(.95,.965,f'{i}/{len(pages)}',fontsize=9,ha='right',va='top',color='#4D5961')
            add_legend(fig,.925)
            fig.text(.045,.045,'Non-exclusive annotation categories; denominator = all representative genes. No enrichment test.',fontsize=8,color='#4D5961')
            pdf.savefig(fig,facecolor='white')
            for ext in ('png','svg'):
                dest=path.with_name(path.name+f'_page{i}').with_suffix('.'+ext)
                fig.savefig(dest,dpi=min(dpi,240),facecolor='white');files.append(str(dest.resolve()))
            plt.close(fig)
    files.append(str(path.with_suffix('.pdf').resolve()))
    return files


def write_selected(path, records, totals):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["plot_order", "category_id", "category_name", "group", "species", "gene_count", "total_genes", "percent_total_genes"])
        for i, row in enumerate(records, 1):
            for sp in SPECIES:
                writer.writerow([i, ";".join(row["key"]), row["label"], row["group"], sp,
                                 row["counts"][sp], totals[sp], f"{row['values'][sp]:.6f}"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tables", type=Path, required=True, help="Directory containing the four verified annotation TSVs")
    parser.add_argument("--out", type=Path, required=True, help="Output directory for PNG, SVG, PDF, selected data and figure metadata")
    parser.add_argument("--go-top", type=int, default=5, help="Main-figure categories per GO ontology (default 5)")
    parser.add_argument("--kegg-top", type=int, default=14, help="Main-figure KEGG level-2 categories (default 14)")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--skip-supplements", action="store_true")
    args = parser.parse_args(argv)
    if not 1 <= args.go_top <= 6 or not 1 <= args.kegg_top <= 15 or args.dpi < 72:
        parser.error("Use --go-top 1–6, --kegg-top 1–15 and --dpi >=72 for this figure layout.")
    try:
        data = load_data(args.tables)
        go = category_records(data["tables"]["goslim_counts.tsv"], ("GO_ID",), "GO_term", "ontology")
        kegg = category_records(data["tables"]["kegg_level2_counts.tsv"], ("A_id", "B_id"), "B_name", "A_name")
        rank = lambda r: (-r["mean"], r["label"], r["key"])
        go_all = [r for ontology in ONTOLOGIES for r in sorted(go, key=rank) if r["group"] == ontology]
        go_main = [r for ontology in ONTOLOGIES for r in sorted(
            [r for r in go if r["group"] == ontology and r["mean"] > 0], key=rank)[:args.go_top]]
        kegg_main = sorted([r for r in kegg if r["mean"] > 0 and r["key"][0] != "09160"
                            and "human disease" not in r["group"].lower()], key=rank)[:args.kegg_top]
        kegg_all = sorted(kegg, key=lambda r: (r["key"][0], -r["mean"], r["label"]))
        if not go_main or not kegg_main:
            raise DataError("No nonzero GO slim or eligible KEGG category data; no main figure can be drawn.")
        plt, font = configure_matplotlib()
        args.out.mkdir(parents=True, exist_ok=True)
        files = main_figure(data, go_main, kegg_main, args.out, plt, args.dpi)
        if not args.skip_supplements:
            files += supplement(go_all, "Complete GO slim annotation profile", args.out / "FigureS_GOslim_all", plt, args.dpi, ONTOLOGY_NAMES)
            files += supplement(kegg_all, "Complete KEGG level-2 annotation profile", args.out / "FigureS_KEGG_level2_all", plt, args.dpi)
        write_selected(args.out / "main_GOslim_selected.tsv", go_main, data["totals"])
        write_selected(args.out / "main_KEGG_selected.tsv", kegg_main, data["totals"])
        manifest = {
            "input_tables_directory": str(args.tables.resolve()), "totals": data["totals"],
            "source_order": data["sources"], "intersection_display_order": data["patterns"],
            "main_GO_categories_per_ontology_max": args.go_top, "main_KEGG_categories_max": args.kegg_top,
            "category_ranking": "descending unweighted mean percentage of all genes across the two species",
            "main_KEGG_exclusion": "human disease categories (A_id=09160 or A_name contains human disease)",
            "percent_denominator": "all longest representative proteins, one per gene",
            "enrichment_test": False, "GO_and_KEGG_categories_nonexclusive": True,
            "figure_size_inches": [7.2, 10.8], "font_family": font,
            "validation": "table schema, counts, denominators, percentages, source coding and exact-intersection sums checked",
            "visual_QA": "not assessed automatically; inspect the rendered figures before publication",
            "files": files,
        }
        (args.out / "functional_figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"Saved {len(files)} figure files to {args.out.resolve()}")
        return 0
    except (DataError, ImportError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
