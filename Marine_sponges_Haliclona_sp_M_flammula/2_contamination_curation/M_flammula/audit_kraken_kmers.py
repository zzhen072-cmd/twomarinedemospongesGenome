#!/usr/bin/env python3
"""Summarise Kraken2 per-k-mer domain support for a selected contig set."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DOMAIN_TAXIDS = {
    2: "Bacteria",
    2157: "Archaea",
    2759: "Eukaryota",
    10239: "Viruses",
}


def canonical(value: str) -> str:
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def load_parent(path: Path) -> dict[int, int]:
    parent: dict[int, int] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split("|")
            if len(fields) >= 2:
                parent[int(fields[0].strip())] = int(fields[1].strip())
    return parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kraken", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--ids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    wanted = {
        canonical(line)
        for line in args.ids.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    parent = load_parent(args.nodes)
    domain_cache: dict[int, str] = {}

    def domain(taxid: int) -> str:
        original = taxid
        visited: list[int] = []
        while taxid not in domain_cache:
            visited.append(taxid)
            if taxid in DOMAIN_TAXIDS:
                domain_cache[taxid] = DOMAIN_TAXIDS[taxid]
                break
            next_taxid = parent.get(taxid, taxid)
            if next_taxid == taxid:
                domain_cache[taxid] = "Other"
                break
            taxid = next_taxid
        result = domain_cache[taxid]
        for item in visited:
            domain_cache[item] = result
        return domain_cache.get(original, result)

    rows: list[dict[str, object]] = []
    found: set[str] = set()
    with args.kraken.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Kraken line {line_number} has fewer than five columns")
            contig_id = canonical(fields[1])
            if contig_id not in wanted:
                continue
            if contig_id in found:
                raise ValueError(f"Duplicate Kraken record for {contig_id}")
            found.add(contig_id)
            counts = {
                "Unclassified": 0,
                "Bacteria": 0,
                "Archaea": 0,
                "Eukaryota": 0,
                "Viruses": 0,
                "Other": 0,
            }
            for token in fields[4].split():
                if ":" not in token:
                    continue
                raw_taxid, raw_count = token.rsplit(":", 1)
                if not raw_taxid.isdigit() or not raw_count.isdigit():
                    continue
                taxid = int(raw_taxid)
                count = int(raw_count)
                if taxid == 0:
                    counts["Unclassified"] += count
                else:
                    counts[domain(taxid)] += count
            classified = sum(value for key, value in counts.items() if key != "Unclassified")
            total = classified + counts["Unclassified"]
            prokaryotic = counts["Bacteria"] + counts["Archaea"]
            rows.append({
                "contig_id": contig_id,
                "kraken_status": fields[0],
                "kraken_call_taxid": fields[2],
                "sequence_length_bp": fields[3],
                "kmer_total": total,
                "kmer_unclassified": counts["Unclassified"],
                "kmer_bacteria": counts["Bacteria"],
                "kmer_archaea": counts["Archaea"],
                "kmer_eukaryota": counts["Eukaryota"],
                "kmer_viruses": counts["Viruses"],
                "kmer_other": counts["Other"],
                "classified_kmer_count": classified,
                "classified_prokaryote_fraction": (
                    f"{prokaryotic / classified:.12g}" if classified else ""
                ),
                "all_kmer_prokaryote_fraction": (
                    f"{prokaryotic / total:.12g}" if total else ""
                ),
                "classified_eukaryote_fraction": (
                    f"{counts['Eukaryota'] / classified:.12g}" if classified else ""
                ),
            })

    missing = wanted - found
    if missing:
        raise ValueError(f"Missing {len(missing)} selected IDs in Kraken output")
    rows.sort(key=lambda row: row["contig_id"])
    fields = list(rows[0]) if rows else []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} records to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
