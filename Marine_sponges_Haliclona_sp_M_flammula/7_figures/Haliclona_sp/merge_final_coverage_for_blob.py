#!/usr/bin/env python3
"""Replace inherited LHT blob-plot depth with audited final-assembly depth."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def index_unique(rows: list[dict[str, str]], source: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, 2):
        ident = row.get("contig_id", "").strip()
        if not ident:
            raise ValueError(f"{source}:{line_number}: blank contig_id")
        if ident in result:
            raise ValueError(f"{source}:{line_number}: duplicate contig_id {ident}")
        result[ident] = row
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()

    base_fields, base_rows = read_tsv(args.base)
    coverage_fields, coverage_rows = read_tsv(args.coverage)
    required_base = {"contig_id", "length_bp", "gc_fraction", "ont_coverage", "phylum"}
    required_coverage = {
        "sample_id",
        "contig_id",
        "length_bp",
        "gc_fraction",
        "default_mean_depth_x",
        "primary_mean_depth_x",
        "primary_q20_mean_depth_x",
    }
    if not required_base.issubset(base_fields):
        raise ValueError(f"Base table lacks columns: {sorted(required_base - set(base_fields))}")
    if not required_coverage.issubset(coverage_fields):
        raise ValueError(
            f"Coverage table lacks columns: {sorted(required_coverage - set(coverage_fields))}"
        )

    base_by_id = index_unique(base_rows, str(args.base))
    coverage_by_id = index_unique(coverage_rows, str(args.coverage))
    if set(base_by_id) != set(coverage_by_id):
        raise ValueError(
            "Contig ID mismatch: "
            f"base-only={len(set(base_by_id) - set(coverage_by_id))}, "
            f"coverage-only={len(set(coverage_by_id) - set(base_by_id))}"
        )

    merged: list[dict[str, str]] = []
    max_gc_difference = 0.0
    zero_contigs = 0
    zero_span = 0
    for row in base_rows:
        ident = row["contig_id"].strip()
        coverage = coverage_by_id[ident]
        if coverage["sample_id"] != "LHT-04":
            raise ValueError(f"Unexpected sample_id for {ident}: {coverage['sample_id']}")
        base_length = int(row["length_bp"])
        coverage_length = int(coverage["length_bp"])
        if base_length != coverage_length:
            raise ValueError(
                f"Length mismatch for {ident}: {base_length} != {coverage_length}"
            )
        gc_difference = abs(float(row["gc_fraction"]) - float(coverage["gc_fraction"]))
        max_gc_difference = max(max_gc_difference, gc_difference)
        if gc_difference > 1e-12:
            raise ValueError(f"GC mismatch for {ident}: difference={gc_difference:.3g}")
        depth = float(coverage["default_mean_depth_x"])
        if not math.isfinite(depth) or depth < 0:
            raise ValueError(f"Invalid default depth for {ident}: {depth}")
        if depth == 0:
            zero_contigs += 1
            zero_span += base_length
        updated = dict(row)
        updated["ont_coverage"] = f"{depth:.15g}"
        merged.append(updated)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=base_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(merged)

    audit = {
        "status": "PASS",
        "coverage_definition": "default compatibility depth; mosdepth -F 1796 -Q 0",
        "base_table_sha256": sha256(args.base),
        "coverage_table_sha256": sha256(args.coverage),
        "output_table_sha256": sha256(args.output),
        "contigs": len(merged),
        "span_bp": sum(int(row["length_bp"]) for row in merged),
        "id_sets_equal": True,
        "lengths_equal": True,
        "max_gc_fraction_difference": max_gc_difference,
        "fully_zero_coverage_contigs": zero_contigs,
        "fully_zero_coverage_contig_span_bp": zero_span,
        "excluded_rows": 0,
    }
    args.audit.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
