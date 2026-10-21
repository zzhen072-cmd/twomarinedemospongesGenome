#!/usr/bin/env python3
"""Audit BRAKER gene models and projected metazoan BUSCO changes before removal."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import re
from datetime import datetime, timezone
from pathlib import Path


def canonical(value: str) -> str:
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--contig-map", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--busco-table", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    short_to_original: dict[str, str] = {}
    with args.contig_map.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Contig-map line {line_number} has fewer than two columns")
            short, original = fields[:2]
            short_to_original[short] = canonical(original)

    contig_genes: dict[str, set[str]] = collections.defaultdict(set)
    contig_introns: collections.Counter[str] = collections.Counter()
    gene_to_contig: dict[str, str] = {}
    with gzip.open(args.gff, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] not in short_to_original:
                continue
            contig = short_to_original[fields[0]]
            feature = fields[2]
            attributes = fields[8]
            if feature == "gene":
                match = re.search(r"(?:^|;)ID=([^;]+)", attributes)
                if match:
                    gene = match.group(1)
                    contig_genes[contig].add(gene)
                    gene_to_contig[gene] = contig
            elif feature == "intron":
                contig_introns[contig] += 1

    busco_hits: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    all_busco_ids: set[str] = set()
    with args.busco_table.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            busco_id = fields[0]
            status = fields[1]
            all_busco_ids.add(busco_id)
            if status == "Missing" or len(fields) < 3:
                continue
            sequence = fields[2]
            gene = sequence.split(".", 1)[0]
            busco_hits[busco_id].append({
                "status": status,
                "sequence": sequence,
                "gene": gene,
                "contig": gene_to_contig.get(gene, ""),
            })

    with args.review.open(encoding="utf-8-sig", newline="") as handle:
        review = list(csv.DictReader(handle, delimiter="\t"))
    decisions = {canonical(row["contig_id"]): row["final_decision"] for row in review}
    removed = {ident for ident, decision in decisions.items() if decision == "remove_confirmed_contaminant"}

    contig_busco_rows: dict[str, list[tuple[str, str, str]]] = collections.defaultdict(list)
    for busco_id, hits in busco_hits.items():
        for hit in hits:
            if hit["contig"]:
                contig_busco_rows[hit["contig"]].append(
                    (busco_id, hit["status"], hit["sequence"])
                )

    rows: list[dict[str, object]] = []
    for source in review:
        ident = canonical(source["contig_id"])
        if source["final_decision"] not in {
            "remove_confirmed_contaminant", "retain_host_supported",
            "retain_ambiguous_review"
        }:
            continue
        length = int(source["length_bp"])
        gene_count = len(contig_genes.get(ident, set()))
        intron_count = contig_introns[ident]
        hits = sorted(contig_busco_rows.get(ident, []))
        rows.append({
            "contig_id": ident,
            "length_bp": length,
            "final_decision": source["final_decision"],
            "evidence_class": source["evidence_class"],
            "braker_gene_count": gene_count,
            "braker_genes_per_mb": f"{gene_count / length * 1_000_000:.12g}",
            "braker_intron_count": intron_count,
            "braker_introns_per_gene": f"{intron_count / gene_count:.12g}" if gene_count else "",
            "metazoa_busco_hit_count": len(hits),
            "metazoa_busco_ids": ";".join(hit[0] for hit in hits),
            "metazoa_busco_statuses": ";".join(hit[1] for hit in hits),
            "metazoa_busco_sequences": ";".join(hit[2] for hit in hits),
        })
    fields = list(rows[0])
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    def projected_status(busco_id: str, remove_set: set[str]) -> str:
        hits = [hit for hit in busco_hits.get(busco_id, []) if hit["contig"] not in remove_set]
        complete_hits = [hit for hit in hits if hit["status"] in {"Complete", "Duplicated"}]
        if len(complete_hits) >= 2:
            return "Duplicated"
        if len(complete_hits) == 1:
            return "Complete"
        if any(hit["status"] == "Fragmented" for hit in hits):
            return "Fragmented"
        return "Missing"

    original_status = {busco_id: projected_status(busco_id, set()) for busco_id in all_busco_ids}
    filtered_status = {busco_id: projected_status(busco_id, removed) for busco_id in all_busco_ids}
    original_counts = collections.Counter(original_status.values())
    filtered_counts = collections.Counter(filtered_status.values())
    changed = [
        {
            "busco_id": busco_id,
            "before": original_status[busco_id],
            "projected_after": filtered_status[busco_id],
            "removed_hits": [
                hit for hit in busco_hits.get(busco_id, []) if hit["contig"] in removed
            ],
            "retained_hits": [
                hit for hit in busco_hits.get(busco_id, []) if hit["contig"] not in removed
            ],
        }
        for busco_id in sorted(all_busco_ids)
        if original_status[busco_id] != filtered_status[busco_id]
    ]
    total = len(all_busco_ids)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "Projected protein-mode metazoa BUSCO effect using the existing BRAKER "
            "annotation; this is a removal safety audit, not a replacement for rerunning BUSCO."
        ),
        "busco_total": total,
        "original_projected_counts_from_full_table": dict(original_counts),
        "filtered_projected_counts_from_full_table": dict(filtered_counts),
        "original_complete_percent": (
            original_counts["Complete"] + original_counts["Duplicated"]
        ) / total * 100,
        "filtered_complete_percent": (
            filtered_counts["Complete"] + filtered_counts["Duplicated"]
        ) / total * 100,
        "changed_busco_count": len(changed),
        "changed_buscos": changed,
        "removed_contigs_with_any_metazoa_busco_hit": sum(
            row["final_decision"] == "remove_confirmed_contaminant"
            and int(row["metazoa_busco_hit_count"]) > 0
            for row in rows
        ),
        "interpretation": (
            "A small loss can reflect removal of bacterial homologues that spuriously "
            "matched conserved metazoan profiles; rerun genome/protein BUSCO on the "
            "filtered derivative before publication."
        ),
    }
    args.output_json.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in audit.items() if key != "changed_buscos"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
