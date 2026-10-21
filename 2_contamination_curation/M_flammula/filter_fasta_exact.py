#!/usr/bin/env python3
"""Create retained/removed FASTA derivatives with exact-ID conservation checks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
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
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--remove-ids", type=Path, required=True)
    parser.add_argument("--retained-fasta", type=Path, required=True)
    parser.add_argument("--removed-fasta", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    parser.add_argument("--expected-input-count", type=int, required=True)
    parser.add_argument("--expected-input-bp", type=int, required=True)
    parser.add_argument("--expected-remove-count", type=int, required=True)
    parser.add_argument("--expected-remove-bp", type=int, required=True)
    args = parser.parse_args()

    destinations = [args.retained_fasta, args.removed_fasta, args.audit_json]
    existing = [str(path) for path in destinations if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing output: " + ", ".join(existing))
    for path in destinations:
        path.parent.mkdir(parents=True, exist_ok=True)

    remove_ids = {
        canonical(line)
        for line in args.remove_ids.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if len(remove_ids) != args.expected_remove_count:
        raise ValueError(
            f"Remove-ID count is {len(remove_ids)}, expected {args.expected_remove_count}"
        )

    token = f".tmp.{os.getpid()}"
    retained_tmp = args.retained_fasta.with_name(args.retained_fasta.name + token)
    removed_tmp = args.removed_fasta.with_name(args.removed_fasta.name + token)
    seen: set[str] = set()
    found_removed: set[str] = set()
    counts = {"retained": 0, "removed": 0}
    spans = {"retained": 0, "removed": 0}

    current_id: str | None = None
    current_header = ""
    sequence_lines: list[str] = []

    def write_record(retained_handle, removed_handle) -> None:
        nonlocal current_id, current_header, sequence_lines
        if current_id is None:
            return
        if current_id in seen:
            raise ValueError(f"Duplicate FASTA identifier: {current_id}")
        seen.add(current_id)
        sequence_bp = sum(len(line.strip()) for line in sequence_lines)
        key = "removed" if current_id in remove_ids else "retained"
        if key == "removed":
            found_removed.add(current_id)
        handle = removed_handle if key == "removed" else retained_handle
        handle.write(current_header)
        for line in sequence_lines:
            handle.write(line)
        counts[key] += 1
        spans[key] += sequence_bp

    try:
        with args.input_fasta.open(encoding="utf-8") as source, \
                retained_tmp.open("w", encoding="utf-8") as retained_handle, \
                removed_tmp.open("w", encoding="utf-8") as removed_handle:
            for line_number, line in enumerate(source, 1):
                if line.startswith(">"):
                    write_record(retained_handle, removed_handle)
                    current_header = line
                    current_id = canonical(line[1:])
                    sequence_lines = []
                else:
                    if current_id is None:
                        raise ValueError(f"Sequence before first FASTA header at line {line_number}")
                    sequence_lines.append(line)
            write_record(retained_handle, removed_handle)

        missing = remove_ids - found_removed
        unexpected_input_count = counts["retained"] + counts["removed"]
        unexpected_input_bp = spans["retained"] + spans["removed"]
        errors = []
        if missing:
            errors.append(f"{len(missing)} remove IDs were absent from the FASTA")
        if unexpected_input_count != args.expected_input_count:
            errors.append(
                f"input count {unexpected_input_count} != {args.expected_input_count}"
            )
        if unexpected_input_bp != args.expected_input_bp:
            errors.append(f"input span {unexpected_input_bp} != {args.expected_input_bp}")
        if counts["removed"] != args.expected_remove_count:
            errors.append(f"removed count {counts['removed']} != {args.expected_remove_count}")
        if spans["removed"] != args.expected_remove_bp:
            errors.append(f"removed span {spans['removed']} != {args.expected_remove_bp}")
        if errors:
            raise RuntimeError("; ".join(errors))

        os.replace(retained_tmp, args.retained_fasta)
        os.replace(removed_tmp, args.removed_fasta)
    except Exception:
        for temporary in (retained_tmp, removed_tmp):
            if temporary.exists():
                temporary.unlink()
        raise

    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "path": str(args.input_fasta),
            "sha256": sha256(args.input_fasta),
            "contig_count": counts["retained"] + counts["removed"],
            "assembly_span_bp": spans["retained"] + spans["removed"],
        },
        "retained": {
            "path": str(args.retained_fasta),
            "sha256": sha256(args.retained_fasta),
            "contig_count": counts["retained"],
            "assembly_span_bp": spans["retained"],
        },
        "removed": {
            "path": str(args.removed_fasta),
            "sha256": sha256(args.removed_fasta),
            "contig_count": counts["removed"],
            "assembly_span_bp": spans["removed"],
        },
        "checks": {
            "remove_ids_unique": True,
            "all_remove_ids_found_once": True,
            "retained_removed_disjoint": True,
            "count_conservation": True,
            "span_conservation": True,
            "input_was_not_modified": True,
        },
    }
    args.audit_json.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
