#!/usr/bin/env python3
"""Build an audited LAN GC/ONT-coverage table with BLAST-derived phyla."""

from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path


def canonical_id(value: str) -> str:
    """Use the BLAST-style lcl| prefix while accepting either input form."""
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows, fields=None) -> None:
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ValueError(f"Cannot infer fields for empty table: {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def fasta_stats(path: Path) -> dict[str, dict]:
    result: dict[str, dict] = {}
    ident = None
    length = gc = acgt = n_bases = 0

    def save() -> None:
        if ident is None:
            return
        if ident in result:
            raise ValueError(f"Duplicate FASTA identifier: {ident}")
        result[ident] = {
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
            else:
                sequence = line.strip().upper()
                length += len(sequence)
                gc += sequence.count("G") + sequence.count("C")
                acgt += sum(sequence.count(base) for base in "ACGT")
                n_bases += sequence.count("N")
    save()
    return result


def read_coverage_table(path: Path) -> dict[str, dict]:
    required = {
        "contig_id",
        "length_bp",
        "gc_fraction",
        "n_bases",
        "ONT_mean_depth",
    }
    result = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Coverage table lacks columns: " + ", ".join(sorted(missing)))
        for line_number, row in enumerate(reader, 2):
            ident = canonical_id(row["contig_id"])
            if not ident or ident in result:
                raise ValueError(f"Line {line_number}: empty or duplicate contig_id")
            length = int(row["length_bp"])
            gc_fraction = float(row["gc_fraction"])
            coverage = float(row["ONT_mean_depth"])
            if (
                length <= 0
                or not 0 <= gc_fraction <= 1
                or not math.isfinite(coverage)
                or coverage < 0
            ):
                raise ValueError(f"Line {line_number}: invalid length, GC, or ONT depth")
            result[ident] = {
                "length_bp": length,
                "gc_fraction": gc_fraction,
                "n_bases": int(row["n_bases"]),
                "ont_coverage": coverage,
            }
    return result


def union_length(intervals) -> int:
    total = 0
    end = 0
    for low, high in sorted(intervals):
        start = max(low, end + 1)
        if high >= start:
            total += high - start + 1
        end = max(end, high)
    return total


def n50(lengths) -> int:
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def summarize(rows):
    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[row["phylum"]].append(row["length_bp"])
    total = sum(row["length_bp"] for row in rows)
    return [
        {
            "phylum": label,
            "contig_count": len(lengths),
            "sum_length_bp": sum(lengths),
            "length_fraction": sum(lengths) / total,
            "n50_bp": n50(lengths),
        }
        for label, lengths in sorted(grouped.items(), key=lambda item: -sum(item[1]))
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument(
        "--fasta",
        type=Path,
        help=(
            "Optional final FASTA for an additional independent sequence check. "
            "When omitted, sequence metrics come from the audited coverage table."
        ),
    )
    parser.add_argument("--coverage-table", type=Path, required=True)
    parser.add_argument("--taxonomy-json", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--blast-complete",
        action="store_true",
        help="Confirm that the BLAST file represents the completed whole-assembly search.",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    coverage = read_coverage_table(args.coverage_table)
    if args.fasta:
        assembly = fasta_stats(args.fasta)
        if set(assembly) != set(coverage):
            raise ValueError(
                f"FASTA/coverage ID mismatch: FASTA-only={len(set(assembly)-set(coverage))}, "
                f"coverage-only={len(set(coverage)-set(assembly))}"
            )
        sequence_mismatches = []
        for ident in assembly:
            if (
                assembly[ident]["length_bp"] != coverage[ident]["length_bp"]
                or assembly[ident]["n_bases"] != coverage[ident]["n_bases"]
                or abs(assembly[ident]["gc_fraction"] - coverage[ident]["gc_fraction"])
                > 1e-12
            ):
                sequence_mismatches.append(ident)
        if sequence_mismatches:
            raise ValueError(
                "Coverage table length/GC/N metrics differ from FASTA for "
                f"{len(sequence_mismatches)} contigs"
            )
    else:
        assembly = {
            ident: {
                "length_bp": row["length_bp"],
                "gc_fraction": row["gc_fraction"],
                "n_bases": row["n_bases"],
            }
            for ident, row in coverage.items()
        }

    taxonomy = json.loads(args.taxonomy_json.read_text(encoding="utf-8"))
    raw_rows = list(csv.reader(args.blast.open(encoding="utf-8"), delimiter="\t"))
    hits = collections.defaultdict(list)
    annotated = []
    duplicates = []
    seen = {}
    for line_number, row in enumerate(raw_rows, 1):
        if len(row) not in (14, 15):
            raise ValueError(f"Unexpected BLAST column count at line {line_number}: {len(row)}")
        ident = canonical_id(row[0])
        if canonical_id(row[3]) != ident:
            raise ValueError(f"BLAST repeated query ID differs at line {line_number}")
        if ident not in assembly:
            raise ValueError(f"BLAST query absent from final assembly: {ident}")
        low, high = sorted((int(row[9]), int(row[10])))
        if not 1 <= low <= high <= assembly[ident]["length_bp"]:
            raise ValueError(f"BLAST query coordinates out of range at line {line_number}")
        requested_taxids = row[1].split(";")
        try:
            taxon_records = [taxonomy[taxid] for taxid in requested_taxids]
        except KeyError as exc:
            raise ValueError(f"Missing cached taxonomy for TaxID {exc.args[0]}") from exc
        phyla = {record["phylum"] for record in taxon_records}
        if len(phyla) == 1 and "" not in phyla:
            phylum = next(iter(phyla))
            resolution = "same_phylum" if len(requested_taxids) > 1 else "single_taxid"
        else:
            phylum = "unresolved"
            resolution = "multiple_phyla_or_no_phylum"
        key = tuple(row)
        duplicate = key in seen
        if duplicate:
            duplicates.append(
                {
                    "duplicate_line": line_number,
                    "first_line": seen[key],
                    "contig_id": ident,
                }
            )
        else:
            seen[key] = line_number
        hit = {
            "contig_id": ident,
            "line_number": line_number,
            "taxids": row[1],
            "subject_accession": row[4],
            "bitscore": float(row[2]),
            "pident": float(row[5]),
            "alignment_length": int(row[6]),
            "qstart": low,
            "qend": high,
            "evalue": float(row[13]),
            "phylum": phylum,
            "taxid_resolution": resolution,
            "duplicate": duplicate,
        }
        hits[ident].append(hit)
        annotated.append(hit)

    merged = []
    for ident, metrics in coverage.items():
        raw_hits = hits.get(ident, [])
        unique_hits = [hit for hit in raw_hits if not hit["duplicate"]]
        scores = collections.defaultdict(float)
        for hit in unique_hits:
            scores[hit["phylum"]] += hit["bitscore"]
        if scores:
            top_score = max(scores.values())
            winners = sorted(
                label
                for label, score in scores.items()
                if math.isclose(score, top_score, rel_tol=1e-12)
            )
            phylum = winners[0] if len(winners) == 1 else "unresolved"
            status = (
                "assigned"
                if phylum != "unresolved"
                else ("tied_phylum_scores" if len(winners) > 1 else "unresolved_phylum")
            )
            best = sorted(unique_hits, key=lambda hit: (-hit["bitscore"], hit["line_number"]))[0]
        else:
            phylum = "no-hit" if args.blast_complete else "not-in-supplied-BLAST"
            status = (
                "no_hit_in_completed_search"
                if args.blast_complete
                else "search_completeness_unconfirmed"
            )
            top_score = 0.0
            best = {}
        total_score = sum(scores.values())
        merged.append(
            {
                "contig_id": ident,
                **metrics,
                "n_bases": assembly[ident]["n_bases"],
                "phylum": phylum,
                "assignment_status": status,
                "raw_blast_rows": len(raw_hits),
                "unique_blast_rows": len(unique_hits),
                "duplicate_rows_removed": len(raw_hits) - len(unique_hits),
                "phylum_bitscore_sum": top_score if scores else "",
                "phylum_score_fraction": top_score / total_score if total_score else "",
                "phylum_count": len(scores),
                "phylum_scores_json": json.dumps(
                    dict(sorted(scores.items())), ensure_ascii=False, separators=(",", ":")
                ),
                "query_aligned_bp": union_length(
                    (hit["qstart"], hit["qend"]) for hit in unique_hits
                ),
                "best_hit_accession": best.get("subject_accession", ""),
                "best_hit_taxids": best.get("taxids", ""),
                "best_hit_bitscore": best.get("bitscore", ""),
                "best_hit_pident": best.get("pident", ""),
                "best_hit_evalue": best.get("evalue", ""),
                "best_hit_phylum": best.get("phylum", ""),
            }
        )

    write_tsv(args.outdir / "lan_final_plot_data.tsv", merged)
    write_tsv(args.outdir / "lan_blast_hits_annotated.tsv", annotated)
    write_tsv(
        args.outdir / "duplicate_blast_rows.tsv",
        duplicates,
        ["duplicate_line", "first_line", "contig_id"],
    )
    write_tsv(args.outdir / "lan_final_phylum_summary.tsv", summarize(merged))

    total_length = sum(row["length_bp"] for row in merged)
    weighted_depth = sum(
        row["length_bp"] * row["ont_coverage"] for row in merged
    ) / total_length
    audit = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sample": "LAN",
        "assembly_contigs": len(merged),
        "assembly_span_bp": total_length,
        "assembly_n50_bp": n50(row["length_bp"] for row in merged),
        "length_weighted_ont_depth": weighted_depth,
        "zero_coverage_contigs": sum(row["ont_coverage"] == 0 for row in merged),
        "blast_complete_confirmed": args.blast_complete,
        "blast_rows": len(raw_rows),
        "exact_duplicate_rows_removed": len(duplicates),
        "blast_hit_queries": len(hits),
        "no_hit_contigs": sum(row["phylum"] == "no-hit" for row in merged),
        "no_hit_span_bp": sum(
            row["length_bp"] for row in merged if row["phylum"] == "no-hit"
        ),
        "taxonomy_record_count": len(taxonomy),
        "method": (
            "Remove exact duplicate BLAST rows. Resolve multi-TaxID hits to a phylum only "
            "when all TaxIDs share that phylum. Sum bit scores by phylum within each contig; "
            "ties are unresolved. No additional e-value or identity cutoff is applied."
        ),
        "inputs": [
            {
                "role": role,
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for role, path in (
                ("whole-assembly BLASTN table", args.blast),
                ("final-assembly ONT coverage table", args.coverage_table),
                ("cached NCBI Taxonomy records", args.taxonomy_json),
            )
        ]
        + (
            [
                {
                    "role": "final assembly FASTA",
                    "sha256": sha256(args.fasta),
                    "size_bytes": args.fasta.stat().st_size,
                }
            ]
            if args.fasta
            else []
        ),
        "phylum_summary": summarize(merged),
    }
    (args.outdir / "analysis_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in audit.items() if key != "inputs"}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
