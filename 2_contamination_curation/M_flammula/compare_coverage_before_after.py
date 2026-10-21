#!/usr/bin/env python3
"""Compare retained-contig ONT depths before and after reference filtering."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path


def canonical(value: str) -> str:
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def read_table(path: Path) -> dict[str, dict[str, str]]:
    result = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ident = canonical(row["contig_id"])
            if ident in result:
                raise ValueError(f"Duplicate contig ID in {path}: {ident}")
            result[ident] = row
    return result


def percentile(values: list[float], probability: float) -> float:
    values = sorted(values)
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    before = read_table(args.before)
    after = read_table(args.after)
    if not set(after) <= set(before):
        raise ValueError("After table contains IDs absent from the before table")

    rows = []
    old_depths = []
    new_depths = []
    for ident in sorted(after):
        old = before[ident]
        new = after[ident]
        if (
            int(old["length_bp"]) != int(new["length_bp"])
            or int(old["n_bases"]) != int(new["n_bases"])
            or abs(float(old["gc_fraction"]) - float(new["gc_fraction"])) > 1e-12
        ):
            raise ValueError(f"Sequence metrics changed for retained contig {ident}")
        old_depth = float(old["ONT_mean_depth"])
        new_depth = float(new["ONT_mean_depth"])
        difference = new_depth - old_depth
        relative = difference / old_depth if old_depth else None
        old_depths.append(old_depth)
        new_depths.append(new_depth)
        rows.append({
            "contig_id": ident,
            "length_bp": old["length_bp"],
            "gc_fraction": old["gc_fraction"],
            "before_ONT_mean_depth": f"{old_depth:.12g}",
            "after_ONT_mean_depth": f"{new_depth:.12g}",
            "depth_difference": f"{difference:.12g}",
            "relative_depth_change": "" if relative is None else f"{relative:.12g}",
        })

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    differences = [new - old for old, new in zip(old_depths, new_depths)]
    absolute = [abs(value) for value in differences]
    relative_nonzero = [
        abs(new - old) / old
        for old, new in zip(old_depths, new_depths)
        if old > 0
    ]
    mean_old = statistics.fmean(old_depths)
    mean_new = statistics.fmean(new_depths)
    covariance = statistics.fmean(
        (old - mean_old) * (new - mean_new)
        for old, new in zip(old_depths, new_depths)
    )
    variance_old = statistics.fmean((value - mean_old) ** 2 for value in old_depths)
    variance_new = statistics.fmean((value - mean_new) ** 2 for value in new_depths)
    correlation = covariance / math.sqrt(variance_old * variance_new)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "retained_contig_count": len(rows),
        "removed_reference_contig_count": len(before) - len(after),
        "sequence_metrics_identical_for_all_retained_contigs": True,
        "pearson_depth_correlation": correlation,
        "absolute_depth_change": {
            "median": statistics.median(absolute),
            "p95": percentile(absolute, 0.95),
            "p99": percentile(absolute, 0.99),
            "maximum": max(absolute),
        },
        "absolute_relative_depth_change_for_nonzero_before": {
            "median": statistics.median(relative_nonzero),
            "p95": percentile(relative_nonzero, 0.95),
            "p99": percentile(relative_nonzero, 0.99),
            "maximum": max(relative_nonzero),
            "count_over_1_percent": sum(value > 0.01 for value in relative_nonzero),
            "count_over_5_percent": sum(value > 0.05 for value in relative_nonzero),
        },
        "zero_depth": {
            "before": sum(value == 0 for value in old_depths),
            "after": sum(value == 0 for value in new_depths),
        },
    }
    args.output_json.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
