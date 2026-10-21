#!/usr/bin/env python3
"""Filter the audited LAN BLAST table to the retained derivative-assembly IDs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def canonical(value: str) -> str:
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--retained-ids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.audit.exists():
        raise FileExistsError("Refusing to overwrite filtered BLAST output or audit")
    retained = {
        canonical(line)
        for line in args.retained_ids.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    input_rows = output_rows = 0
    input_queries: set[str] = set()
    output_queries: set[str] = set()
    with args.input.open(encoding="utf-8", newline="") as source, \
            args.output.open("w", encoding="utf-8", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        for line_number, row in enumerate(reader, 1):
            if len(row) not in {14, 15}:
                raise ValueError(f"Unexpected BLAST column count at line {line_number}: {len(row)}")
            query = canonical(row[0])
            if canonical(row[3]) != query:
                raise ValueError(f"Repeated query identifier differs at line {line_number}")
            input_rows += 1
            input_queries.add(query)
            if query in retained:
                writer.writerow(row)
                output_rows += 1
                output_queries.add(query)
    if output_queries - retained:
        raise RuntimeError("Filtered BLAST contains a query absent from retained IDs")
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "path": str(args.input), "sha256": sha256(args.input),
            "row_count": input_rows, "query_count": len(input_queries),
        },
        "retained_id_count": len(retained),
        "output": {
            "path": str(args.output), "sha256": sha256(args.output),
            "row_count": output_rows, "query_count": len(output_queries),
        },
        "removed_blast_rows": input_rows - output_rows,
        "removed_hit_queries": len(input_queries - output_queries),
        "checks": {
            "output_queries_subset_of_retained_ids": True,
            "raw_rows_preserved_without_reformatting_except_line_endings": True,
        },
    }
    args.audit.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
