#!/usr/bin/env python3
"""Combine LAN FASTA sequence metrics with whole-contig mosdepth output."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


def canonical_id(value: str) -> str:
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def read_fasta(path: Path) -> dict[str, dict[str, int | float]]:
    records: dict[str, dict[str, int | float]] = {}
    ident = None
    length = gc = acgt = n_bases = 0

    def save() -> None:
        if ident is None:
            return
        if ident in records:
            raise ValueError(f"Duplicate FASTA identifier: {ident}")
        records[ident] = {
            "length_bp": length,
            "gc_fraction": gc / acgt if acgt else 0.0,
            "n_bases": n_bases,
        }

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                save()
                ident = canonical_id(line[1:])
                length = gc = acgt = n_bases = 0
                continue
            sequence = line.strip().upper()
            length += len(sequence)
            gc += sequence.count("G") + sequence.count("C")
            acgt += sum(sequence.count(base) for base in "ACGT")
            n_bases += sequence.count("N")
    save()
    if not records:
        raise ValueError("FASTA contains no sequences")
    return records


def read_regions(path: Path) -> dict[str, dict[str, int | float]]:
    result: dict[str, dict[str, int | float]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) < 4:
                raise ValueError(f"Line {line_number}: expected four mosdepth columns")
            ident = canonical_id(row[0])
            start, end, depth = int(row[1]), int(row[2]), float(row[3])
            if ident in result:
                raise ValueError(f"Line {line_number}: duplicate region for {ident}")
            if start != 0 or end <= 0 or not math.isfinite(depth) or depth < 0:
                raise ValueError(f"Line {line_number}: invalid whole-contig interval or depth")
            result[ident] = {"length_bp": end, "ONT_mean_depth": depth}
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()

    assembly = read_fasta(args.fasta)
    regions = read_regions(args.regions)
    if set(assembly) != set(regions):
        raise ValueError(
            f"FASTA/mosdepth ID mismatch: FASTA-only={len(set(assembly)-set(regions))}, "
            f"mosdepth-only={len(set(regions)-set(assembly))}"
        )
    length_mismatches = [
        ident
        for ident in assembly
        if assembly[ident]["length_bp"] != regions[ident]["length_bp"]
    ]
    if length_mismatches:
        raise ValueError(f"Length mismatch for {len(length_mismatches)} contigs")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "contig_id",
                "length_bp",
                "gc_fraction",
                "n_bases",
                "ONT_mean_depth",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for ident, metrics in assembly.items():
            writer.writerow(
                {
                    "contig_id": ident,
                    **metrics,
                    "ONT_mean_depth": regions[ident]["ONT_mean_depth"],
                }
            )

    total_length = sum(int(row["length_bp"]) for row in assembly.values())
    weighted_depth = sum(
        int(assembly[ident]["length_bp"]) * float(regions[ident]["ONT_mean_depth"])
        for ident in assembly
    ) / total_length
    audit = {
        "contig_count": len(assembly),
        "assembly_span_bp": total_length,
        "zero_depth_contigs": sum(
            float(regions[ident]["ONT_mean_depth"]) == 0 for ident in assembly
        ),
        "length_weighted_ONT_depth": weighted_depth,
        "method": "mosdepth -n --by one full-length BED interval per BAM reference",
    }
    args.audit.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
