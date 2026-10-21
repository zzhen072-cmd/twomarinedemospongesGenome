#!/usr/bin/env python3
"""Build figure-specific functional-annotation tables from audited gene rows.

The manuscript reports the InterProScan *any member-database signature*
criterion.  This script therefore derives exact three-way intersections from
the per-gene audit table instead of reusing the narrower integrated-InterPro
intersection table produced for the general annotation delivery.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


SAMPLES = ("LHT-04", "NJ-lan")
SOURCES = (
    ("eggNOG", "eggNOG"),
    ("InterPro", "InterProScan_any_signature"),
    ("SwissProt", "SwissProt"),
)
EXPECTED_TOTALS = {"LHT-04": 23379, "NJ-lan": 35461}
EXPECTED_UNION = {"LHT-04": 19210, "NJ-lan": 26863}


def as_bool(value: str) -> bool:
    value = value.strip().lower()
    if value not in {"true", "false"}:
        raise ValueError(f"Unexpected Boolean value: {value!r}")
    return value == "true"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy_table(source: Path, destination: Path) -> None:
    destination.write_bytes(source.read_bytes())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    gene_table = args.source / "gene_functional_annotations.tsv"
    args.out.mkdir(parents=True, exist_ok=True)

    counts = {sample: Counter() for sample in SAMPLES}
    totals = Counter()
    with gene_table.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species", *(column for _, column in SOURCES)}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing gene-table columns: {sorted(missing)}")
        for row in reader:
            sample = row["species"]
            if sample not in counts:
                raise ValueError(f"Unexpected sample: {sample!r}")
            pattern = "".join("1" if as_bool(row[column]) else "0" for _, column in SOURCES)
            counts[sample][pattern] += 1
            totals[sample] += 1

    if dict(totals) != EXPECTED_TOTALS:
        raise ValueError(f"Unexpected representative-gene totals: {dict(totals)}")

    patterns = [f"{value:03b}" for value in range(8)]
    intersections: list[dict] = []
    coverage: list[dict] = []
    set_sizes: list[dict] = []
    for sample in SAMPLES:
        total = totals[sample]
        union = total - counts[sample]["000"]
        if union != EXPECTED_UNION[sample]:
            raise ValueError(f"{sample}: union {union} != expected {EXPECTED_UNION[sample]}")
        for pattern in patterns:
            active = [name for bit, (name, _) in zip(pattern, SOURCES) if bit == "1"]
            count = counts[sample][pattern]
            intersections.append({
                "species": sample,
                "pattern": pattern,
                "sources": ";".join(active) if active else "None",
                "source_order": ";".join(name for name, _ in SOURCES),
                "gene_count": count,
                "total_genes": total,
                "percent_total_genes": f"{100 * count / total:.6f}",
            })
        for index, (name, _) in enumerate(SOURCES):
            count = sum(counts[sample][pattern] for pattern in patterns if pattern[index] == "1")
            row = {
                "species": sample,
                "source": name,
                "gene_count": count,
                "total_genes": total,
                "percent_total_genes": f"{100 * count / total:.6f}",
            }
            coverage.append(row)
            set_sizes.append(dict(row))

    coverage_fields = ["species", "source", "gene_count", "total_genes", "percent_total_genes"]
    write_tsv(args.out / "functional_coverage.tsv", coverage_fields, coverage)
    write_tsv(args.out / "upset_set_sizes.tsv", coverage_fields, set_sizes)
    write_tsv(
        args.out / "upset_intersections.tsv",
        ["species", "pattern", "sources", "source_order", "gene_count", "total_genes", "percent_total_genes"],
        intersections,
    )
    for name in ("goslim_counts.tsv", "kegg_level2_counts.tsv"):
        copy_table(args.source / name, args.out / name)

    audit = {
        "source_gene_table": str(gene_table.resolve()),
        "samples": list(SAMPLES),
        "totals": dict(totals),
        "source_definition": {
            "eggNOG": "at least one eggNOG-mapper annotation row",
            "InterPro": "at least one InterProScan member-database signature",
            "SwissProt": "at least one retained Swiss-Prot DIAMOND match",
        },
        "exclusive_intersections": {sample: dict(sorted(counts[sample].items())) for sample in SAMPLES},
        "union_counts": {sample: totals[sample] - counts[sample]["000"] for sample in SAMPLES},
        "expected_manuscript_union_counts": EXPECTED_UNION,
        "status": "PASS",
    }
    (args.out / "figure_table_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
