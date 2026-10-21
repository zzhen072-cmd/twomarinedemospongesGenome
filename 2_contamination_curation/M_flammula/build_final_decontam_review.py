#!/usr/bin/env python3
"""Build an audited LAN whole-contig removal decision from independent evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


TIER1 = "tier1_high_confidence_quarantine_candidate"
TIER2 = "tier2_manual_review"
REMOVE = "remove_confirmed_contaminant"
RETAIN_HOST = "retain_host_supported"
RETAIN_REVIEW = "retain_ambiguous_review"
RETAIN = "retain"

PROKARYOTIC_PHYLA = {
    "Acidobacteriota", "Actinomycetota", "Bacteroidota", "Bdellovibrionota",
    "Binatota", "Caldifarciminota", "Chloroflexota", "Cyanobacteriota",
    "Pseudomonadota", "Rhodothermota", "Verrucomicrobiota",
}


def canonical(value: str) -> str:
    value = value.strip().split()[0]
    return value if value.startswith("lcl|") else f"lcl|{value}"


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows, fields


def keyed(rows: list[dict[str, str]], field: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        ident = canonical(row[field])
        if ident in result:
            raise ValueError(f"Duplicate identifier in {field}: {ident}")
        result[ident] = row
    return result


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def n50(lengths: list[int]) -> int:
    target = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict-screen", type=Path, required=True)
    parser.add_argument("--legacy-evidence", type=Path, required=True)
    parser.add_argument("--kraken-kmers", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    screen, screen_fields = read_tsv(args.strict_screen)
    legacy_rows, _ = read_tsv(args.legacy_evidence)
    kraken_rows, _ = read_tsv(args.kraken_kmers)
    legacy = keyed(legacy_rows, "contig")
    kraken = keyed(kraken_rows, "contig_id")
    screen_ids = {canonical(row["contig_id"]) for row in screen}
    if screen_ids - set(legacy):
        raise ValueError(f"Legacy evidence is missing {len(screen_ids - set(legacy))} contigs")
    tier2_ids = {
        canonical(row["contig_id"])
        for row in screen
        if row["strict_tier"] == TIER2
    }
    if set(kraken) != tier2_ids:
        raise ValueError(
            "Kraken audit must contain exactly the Tier 2 IDs: "
            f"missing={len(tier2_ids - set(kraken))}, extra={len(set(kraken) - tier2_ids)}"
        )

    output_rows: list[dict[str, object]] = []
    for source in screen:
        ident = canonical(source["contig_id"])
        old = legacy[ident]
        kk = kraken.get(ident, {})
        length = int(source["length_bp"])
        breadth = float(source["blast_aligned_fraction"] or 0)
        score_fraction = float(source["phylum_score_fraction"] or 0)
        best_identity = float(source["best_hit_pident"] or 0)
        fcs_removed = old["fcs_gx_removed"] == "True"
        tiara_stage1 = old["tiara_stage1"]
        tiara_prokaryote = tiara_stage1 in {"bacteria", "archaea", "prokarya"}
        phylum_scores = json.loads(source["phylum_scores_json"] or "{}")
        resolved_hit_phyla = set(phylum_scores) - {"unresolved"}
        broad_prokaryotic_homology = (
            source["strict_tier"] == TIER2
            and breadth >= 0.50
            and bool(resolved_hit_phyla)
            and resolved_hit_phyla <= PROKARYOTIC_PHYLA
        )
        if source["strict_tier"] == TIER1:
            decision = REMOVE
            evidence_class = "tier1_multiaxis_strict"
            reason = (
                "Pseudomonadota BLAST breadth >=0.50, score fraction >=0.90, "
                "and both GC and ONT depth above Porifera P99"
            )
        elif source["strict_tier"] == TIER2 and fcs_removed:
            decision = REMOVE
            evidence_class = "tier2_blast_plus_fcs_gx"
            reason = (
                "Independent FCS-GX removal flag plus a Pseudomonadota-dominant BLAST "
                "assignment on this short contig"
            )
        elif source["strict_tier"] == TIER2 and tiara_prokaryote:
            decision = REMOVE
            evidence_class = "tier2_blast_plus_tiara_prokaryote"
            reason = (
                "Pseudomonadota BLAST assignment plus independent Tiara whole-sequence "
                f"{tiara_stage1} classification"
            )
        elif broad_prokaryotic_homology:
            decision = REMOVE
            evidence_class = "tier2_broad_prokaryotic_homology"
            reason = (
                "BLAST interval union covers at least 50% of the contig and every "
                "resolved hit phylum is prokaryotic"
            )
        elif source["strict_tier"] == TIER2 and tiara_stage1 == "eukarya":
            decision = RETAIN_HOST
            evidence_class = "tier2_tiara_eukaryotic_support"
            reason = (
                "Tiara supports a eukaryotic whole-sequence origin and the bacterial "
                "BLAST evidence is localized; retain the whole contig"
            )
        elif source["strict_tier"] == TIER2:
            decision = RETAIN_REVIEW
            evidence_class = "tier2_insufficient_independent_evidence"
            reason = (
                "Pseudomonadota label was not corroborated by FCS-GX or Tiara and did "
                "not meet the near-whole-contig homology rule"
            )
        else:
            decision = RETAIN
            evidence_class = "outside_pseudomonadota_screen"
            reason = "No removal decision from the Pseudomonadota-focused review"

        row: dict[str, object] = dict(source)
        row.update({
            "legacy_fcs_gx_removed": fcs_removed,
            "legacy_tiara_stage1": tiara_stage1,
            "legacy_tiara_stage2": old["tiara_stage2"],
            "legacy_decision": old["decision"],
            "kraken_classified_kmer_count": kk.get("classified_kmer_count", ""),
            "kraken_classified_prokaryote_fraction": kk.get(
                "classified_prokaryote_fraction", ""
            ),
            "kraken_classified_eukaryote_fraction": kk.get(
                "classified_eukaryote_fraction", ""
            ),
            "broad_prokaryotic_homology_pass": broad_prokaryotic_homology,
            "final_decision": decision,
            "evidence_class": evidence_class,
            "final_reason": reason,
            "sequence_action": (
                "exclude_from_derivative_host_assembly" if decision == REMOVE
                else "retain_in_derivative_host_assembly"
            ),
        })
        output_rows.append(row)

    extra_fields = [
        "legacy_fcs_gx_removed", "legacy_tiara_stage1", "legacy_tiara_stage2",
        "legacy_decision", "kraken_classified_kmer_count",
        "kraken_classified_prokaryote_fraction",
        "kraken_classified_eukaryote_fraction",
        "broad_prokaryotic_homology_pass", "final_decision", "evidence_class",
        "final_reason", "sequence_action",
    ]
    output_fields = screen_fields + [field for field in extra_fields if field not in screen_fields]
    order = {REMOVE: 0, RETAIN_HOST: 1, RETAIN_REVIEW: 2, RETAIN: 3}
    output_rows.sort(
        key=lambda row: (order[str(row["final_decision"])], -int(row["length_bp"]), row["contig_id"])
    )
    write_tsv(args.outdir / "lan_final_decontam_review.tsv", output_rows, output_fields)

    removed = [row for row in output_rows if row["final_decision"] == REMOVE]
    retained_host = [row for row in output_rows if row["final_decision"] == RETAIN_HOST]
    retained_review = [row for row in output_rows if row["final_decision"] == RETAIN_REVIEW]
    retained = [row for row in output_rows if row["final_decision"] != REMOVE]
    tier2_removed = [row for row in removed if row["strict_tier"] == TIER2]
    for filename, subset in [
        ("remove_confirmed_contaminants.ids.txt", removed),
        ("tier2_removed_confirmed_contaminants.ids.txt", tier2_removed),
        ("tier2_retained_host_supported.ids.txt", retained_host),
        ("tier2_retained_ambiguous_review.ids.txt", retained_review),
        ("retained_derivative_assembly.ids.txt", retained),
    ]:
        (args.outdir / filename).write_text(
            "".join(f"{row['contig_id']}\n" for row in subset), encoding="utf-8"
        )
    write_tsv(args.outdir / "removed_confirmed_contaminants.tsv", removed, output_fields)
    write_tsv(args.outdir / "tier2_removed_confirmed_contaminants.tsv", tier2_removed, output_fields)
    write_tsv(args.outdir / "tier2_retained_host_supported.tsv", retained_host, output_fields)
    write_tsv(args.outdir / "tier2_retained_ambiguous_review.tsv", retained_review, output_fields)

    total_bp = sum(int(row["length_bp"]) for row in output_rows)
    removed_bp = sum(int(row["length_bp"]) for row in removed)
    retained_bp = sum(int(row["length_bp"]) for row in retained)
    by_evidence: dict[str, dict[str, int]] = {}
    for row in removed:
        label = str(row["evidence_class"])
        item = by_evidence.setdefault(label, {"contig_count": 0, "sum_length_bp": 0})
        item["contig_count"] += 1
        item["sum_length_bp"] += int(row["length_bp"])
    summary_rows = [
        {
            "decision": REMOVE,
            "contig_count": len(removed),
            "sum_length_bp": removed_bp,
            "assembly_fraction": f"{removed_bp / total_bp:.12g}",
            "n50_bp": n50([int(row["length_bp"]) for row in removed]),
        },
        {
            "decision": RETAIN_HOST,
            "contig_count": len(retained_host),
            "sum_length_bp": sum(int(row["length_bp"]) for row in retained_host),
            "assembly_fraction": f"{sum(int(row['length_bp']) for row in retained_host) / total_bp:.12g}",
            "n50_bp": n50([int(row["length_bp"]) for row in retained_host]),
        },
        {
            "decision": RETAIN_REVIEW,
            "contig_count": len(retained_review),
            "sum_length_bp": sum(int(row["length_bp"]) for row in retained_review),
            "assembly_fraction": f"{sum(int(row['length_bp']) for row in retained_review) / total_bp:.12g}",
            "n50_bp": n50([int(row["length_bp"]) for row in retained_review]),
        },
        {
            "decision": "retained_derivative_assembly_total",
            "contig_count": len(retained),
            "sum_length_bp": retained_bp,
            "assembly_fraction": f"{retained_bp / total_bp:.12g}",
            "n50_bp": n50([int(row["length_bp"]) for row in retained]),
        },
    ]
    write_tsv(
        args.outdir / "lan_final_decontam_summary.tsv", summary_rows,
        ["decision", "contig_count", "sum_length_bp", "assembly_fraction", "n50_bp"],
    )

    expected = {
        "source": (7366, 302079020),
        "removed": (249, 10189414),
        "tier2_removed": (161, 10071961),
        "tier2_retained_host": (6, 265483),
        "tier2_retained_review": (31, 50282),
        "retained": (7117, 291889606),
    }
    observed = {
        "source": (len(output_rows), total_bp),
        "removed": (len(removed), removed_bp),
        "tier2_removed": (
            len(tier2_removed), sum(int(row["length_bp"]) for row in tier2_removed)
        ),
        "tier2_retained_host": (
            len(retained_host), sum(int(row["length_bp"]) for row in retained_host)
        ),
        "tier2_retained_review": (
            len(retained_review), sum(int(row["length_bp"]) for row in retained_review)
        ),
        "retained": (len(retained), retained_bp),
    }
    if observed != expected:
        raise RuntimeError(f"Review totals changed unexpectedly: {observed!r}")

    large_tier2 = [
        row for row in output_rows
        if row["strict_tier"] == TIER2 and int(row["length_bp"]) >= 100_000
    ]
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "Whole-contig contaminant removal for a sponge host assembly; the original "
            "assembly remains immutable and removed sequences are retained separately."
        ),
        "inputs": {
            "strict_screen": {"path": str(args.strict_screen), "sha256": sha256(args.strict_screen)},
            "legacy_evidence": {"path": str(args.legacy_evidence), "sha256": sha256(args.legacy_evidence)},
            "kraken_kmers": {"path": str(args.kraken_kmers), "sha256": sha256(args.kraken_kmers)},
        },
        "rule": {
            "tier1": "remove all multidimensional strict candidates",
            "tier2_independent": (
                "remove when the Pseudomonadota label is corroborated by a Tiara "
                "prokaryotic whole-sequence call or an independent FCS-GX removal flag"
            ),
            "tier2_broad_prokaryotic_homology": (
                "otherwise remove only when BLAST breadth >=0.50 and every resolved "
                "hit phylum is prokaryotic"
            ),
            "kraken_note": (
                "Kraken classified-k-mer domain fractions are retained as corroborating "
                "audit data but are not a primary removal gate because classified k-mers "
                "are sparse for many records."
            ),
        },
        "expected_and_observed": observed,
        "removed_by_evidence_class": by_evidence,
        "large_tier2_review": {
            "threshold_bp": 100000,
            "contig_count": len(large_tier2),
            "removed_count": sum(row["final_decision"] == REMOVE for row in large_tier2),
            "retained_count": sum(row["final_decision"] != REMOVE for row in large_tier2),
            "retained_ids": [
                row["contig_id"] for row in large_tier2 if row["final_decision"] != REMOVE
            ],
        },
        "safety": {
            "original_assembly_modified": False,
            "removed_sequences_to_be_archived_separately": True,
            "ambiguous_tier2_sequences_removed": False,
        },
    }
    (args.outdir / "decontam_review_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
