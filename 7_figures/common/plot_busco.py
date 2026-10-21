#!/usr/bin/env python3
"""Plot audited BUSCO counts for the locked LHT-04 and NJ-lan assemblies."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/codex-lht-njlan-busco-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from audit_panel_alignment import require_matplotlib_panel_alignment

SAMPLES = ("lht04", "njlan")
DISPLAY = {"lht04": "LHT-04", "njlan": "NJ-lan"}
MODES = ("genome", "proteins")
STATUS = (
    ("single", "Single-copy (S)", "#0072B2"),
    ("duplicated", "Duplicated (D)", "#56B4E9"),
    ("fragmented", "Fragmented (F)", "#E69F00"),
    ("missing", "Missing (M)", "#D7DCE2"),
)

def load_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        raw = list(csv.DictReader(handle, delimiter="\t"))
    selected: list[dict] = []
    for sample in SAMPLES:
        for mode in MODES:
            matches = [row for row in raw if row["sample"] == sample and row["analysis"] == mode]
            if len(matches) != 1:
                raise ValueError(f"Expected one BUSCO row for {sample}/{mode}; found {len(matches)}")
            row = matches[0].copy()
            if row["version"] != "6.0.0" or row["lineage"] != "metazoa_odb10":
                raise ValueError(f"Unexpected BUSCO version or lineage for {sample}/{mode}")
            row["total"] = int(row["total"])
            if row["total"] != 954:
                raise ValueError(f"Unexpected BUSCO denominator for {sample}/{mode}")
            for key, _, _ in STATUS:
                row[key] = int(row[key])
            if sum(row[key] for key, _, _ in STATUS) != row["total"]:
                raise ValueError(f"BUSCO status counts do not sum for {sample}/{mode}")
            row["complete"] = row["single"] + row["duplicated"]
            selected.append(row)
    return selected

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = load_rows(args.summary)

    source_rows = []
    for row in rows:
        output = {
            "sample": DISPLAY[row["sample"]],
            "analysis": "Genome" if row["analysis"] == "genome" else "Representative proteins",
            "BUSCO_version": row["version"], "lineage": row["lineage"],
            "total_groups": row["total"],
        }
        for key, _, _ in STATUS:
            output[f"{key}_count"] = row[key]
            output[f"{key}_percent"] = f"{100 * row[key] / row['total']:.6f}"
        output["complete_count"] = row["complete"]
        output["complete_percent"] = f"{100 * row['complete'] / row['total']:.6f}"
        source_rows.append(output)
    source_path = args.out / "Figure4_BUSCO_source_data.tsv"
    with source_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(source_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(source_rows)

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": .6, "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True)
    fig.subplots_adjust(left=.255, right=.885, bottom=.145, top=.74, hspace=.76)
    fig.text(.04, .962, "BUSCO completeness", fontsize=14, weight="bold", ha="left", va="top")
    fig.text(.04, .90, "BUSCO 6.0.0 | metazoa_odb10 | 954 BUSCO groups", fontsize=9,
             color="#4D5865", ha="left")
    fig.legend(handles=[Patch(facecolor=color, label=label) for _, label, color in STATUS],
               loc="upper left", bbox_to_anchor=(.033, .865), ncol=4, frameon=False,
               handlelength=1.35, handleheight=.85, columnspacing=1.45, fontsize=8.5)

    for ax, sample, letter in zip(axes, SAMPLES, ("a", "b")):
        selected = [row for row in rows if row["sample"] == sample]
        if [row["analysis"] for row in selected] != list(MODES):
            raise ValueError(f"Unexpected BUSCO row order for {sample}")
        for row, y in zip(selected, (.7, 0)):
            left = 0.0
            for key, _, color in STATUS:
                width = 100 * row[key] / row["total"]
                ax.barh(y, width, left=left, height=.40, color=color,
                        edgecolor="white", linewidth=.45, zorder=3)
                left += width
            complete_pct = 100 * row["complete"] / row["total"]
            ax.text(102.5, y, f"{complete_pct:.1f}%", ha="left", va="center",
                    fontsize=9, weight="bold", color="#253449", clip_on=False)
        ax.set_xlim(0, 100); ax.set_ylim(-.4, 1.09)
        ax.set_yticks((.7, 0)); ax.set_yticklabels(("Genome", "Representative\nproteins"))
        ax.tick_params(axis="y", length=0, pad=9, labelsize=9)
        ax.tick_params(axis="x", length=3, width=.6, color="#A7AFB8", labelsize=8)
        ax.grid(axis="x", color="#E6E9ED", linewidth=.6, zorder=0)
        for spine in ("left", "top", "right"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#A7AFB8")
        ax.text(-.34, 1.14, letter, transform=ax.transAxes, fontsize=11, weight="bold")
        ax.text(0, 1.14, DISPLAY[sample], transform=ax.transAxes, fontsize=10, weight="bold")
        ax.text(1.025, 1.14, "C (S + D)", transform=ax.transAxes,
                ha="left", fontsize=8, color="#4D5865")
    axes[-1].set_xticks((0, 20, 40, 60, 80, 100))
    axes[-1].set_xlabel("BUSCO groups (%)", labelpad=6)

    stem = args.out / "Figure4_BUSCO"
    require_matplotlib_panel_alignment(
        fig, axes=axes, panel_ids=("a", "b"), column_groups=(("a", "b"),),
        json_out=stem.with_suffix(".alignment.json"),
        overlay_svg=stem.with_suffix(".alignment.svg"),
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True,
    )
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    manifest = {
        "input_summary": str(args.summary.resolve()), "source_data": str(source_path.resolve()),
        "samples": [DISPLAY[sample] for sample in SAMPLES], "BUSCO_version": "6.0.0",
        "lineage": "metazoa_odb10", "groups": 954,
        "validation": "counts sum to 954; complete recalculated as single + duplicated",
        "protein_mode_independence_caveat": "BRAKER used compleasm/BUSCO hints; protein-mode BUSCO is not fully independent validation",
        "status": "PASS",
    }
    (args.out / "Figure4_BUSCO_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
