#!/usr/bin/env python3
"""Validate the complete LAN evidence-reviewed filtering and blob-plot delivery."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_ids(path: Path) -> set[str]:
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate ID in {path}")
    return set(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    audit_dir = root / "audit"
    assembly_dir = root / "assembly"
    mapping_dir = root / "mapping"
    plot_dir = root / "plot"

    review = read_tsv(audit_dir / "lan_final_decontam_review.tsv")
    if len(review) != 7366 or len({row["contig_id"] for row in review}) != 7366:
        raise RuntimeError("Review table must contain 7,366 unique contigs")
    review_by_id = {row["contig_id"]: row for row in review}
    original_span = sum(int(row["length_bp"]) for row in review)
    removed_ids = read_ids(audit_dir / "remove_confirmed_contaminants.ids.txt")
    retained_ids = read_ids(audit_dir / "retained_derivative_assembly.ids.txt")
    tier2_removed_ids = read_ids(audit_dir / "tier2_removed_confirmed_contaminants.ids.txt")
    tier2_host_ids = read_ids(audit_dir / "tier2_retained_host_supported.ids.txt")
    tier2_retained_ids = read_ids(audit_dir / "tier2_retained_ambiguous_review.ids.txt")
    if removed_ids & retained_ids or removed_ids | retained_ids != set(review_by_id):
        raise RuntimeError("Removed/retained sets are not a disjoint partition of the source")
    removed_span = sum(int(review_by_id[ident]["length_bp"]) for ident in removed_ids)
    retained_span = sum(int(review_by_id[ident]["length_bp"]) for ident in retained_ids)

    filter_audit = json.loads((assembly_dir / "filter_audit.json").read_text(encoding="utf-8"))
    coverage_audit = json.loads((mapping_dir / "coverage_audit.json").read_text(encoding="utf-8"))
    analysis_audit = json.loads((plot_dir / "analysis_audit.json").read_text(encoding="utf-8"))
    blast_audit = json.loads((root / "source_data/blast_filter_audit.json").read_text(encoding="utf-8"))
    plot_rows = read_tsv(plot_dir / "lan_final_plot_data.tsv")
    plot_ids = {row["contig_id"] for row in plot_rows}
    plot_span = sum(int(row["length_bp"]) for row in plot_rows)
    pseudo_rows = [row for row in plot_rows if row["phylum"] == "Pseudomonadota"]

    expected = {
        "source": (7366, 302079020),
        "removed": (249, 10189414),
        "retained": (7117, 291889606),
        "tier2_removed": (161, 10071961),
        "tier2_retained_host": (6, 265483),
        "tier2_retained_review": (31, 50282),
        "retained_pseudomonadota": (37, 315765),
    }
    observed = {
        "source": (len(review), original_span),
        "removed": (len(removed_ids), removed_span),
        "retained": (len(retained_ids), retained_span),
        "tier2_removed": (
            len(tier2_removed_ids),
            sum(int(review_by_id[ident]["length_bp"]) for ident in tier2_removed_ids),
        ),
        "tier2_retained_host": (
            len(tier2_host_ids),
            sum(int(review_by_id[ident]["length_bp"]) for ident in tier2_host_ids),
        ),
        "tier2_retained_review": (
            len(tier2_retained_ids),
            sum(int(review_by_id[ident]["length_bp"]) for ident in tier2_retained_ids),
        ),
        "retained_pseudomonadota": (
            len(pseudo_rows), sum(int(row["length_bp"]) for row in pseudo_rows)
        ),
    }
    if observed != expected:
        raise RuntimeError(f"Unexpected decision or Pseudomonadota totals: {observed!r}")
    if plot_ids != retained_ids or plot_span != retained_span:
        raise RuntimeError("Final plot table is not identical to the retained ID/span set")

    checks = [
        {
            "name": "source_partition_and_expected_totals",
            "passed": True,
            "observed": observed,
        },
        {
            "name": "filter_audit_conservation",
            "passed": (
                filter_audit["input"]["contig_count"] == 7366
                and filter_audit["input"]["assembly_span_bp"] == 302079020
                and filter_audit["removed"]["contig_count"] == 249
                and filter_audit["removed"]["assembly_span_bp"] == 10189414
                and filter_audit["retained"]["contig_count"] == 7117
                and filter_audit["retained"]["assembly_span_bp"] == 291889606
                and all(filter_audit["checks"].values())
            ),
        },
        {
            "name": "remapped_coverage_matches_filtered_assembly",
            "passed": (
                coverage_audit["contig_count"] == 7117
                and coverage_audit["assembly_span_bp"] == 291889606
            ),
            "coverage_audit": coverage_audit,
        },
        {
            "name": "final_analysis_matches_filtered_assembly",
            "passed": (
                analysis_audit["assembly_contigs"] == 7117
                and analysis_audit["assembly_span_bp"] == 291889606
                and analysis_audit["blast_complete_confirmed"] is True
            ),
        },
        {
            "name": "filtered_blast_queries_subset",
            "passed": blast_audit["checks"]["output_queries_subset_of_retained_ids"],
            "output_rows": blast_audit["output"]["row_count"],
            "output_queries": blast_audit["output"]["query_count"],
        },
        {
            "name": "plot_id_set_and_span_match_retained_derivative",
            "passed": plot_ids == retained_ids and plot_span == retained_span,
        },
    ]

    qa_paths = {
        "collision": plot_dir / "lan_review_filtered_blobplot.collision-audit.json",
        "pdf_text": plot_dir / "lan_review_filtered_blobplot.pdf-text-audit.json",
        "alignment": plot_dir / "lan_review_filtered_blobplot.alignment.json",
    }
    if all(path.exists() for path in qa_paths.values()):
        collision = json.loads(qa_paths["collision"].read_text(encoding="utf-8"))
        pdf_text = json.loads(qa_paths["pdf_text"].read_text(encoding="utf-8"))
        alignment = json.loads(qa_paths["alignment"].read_text(encoding="utf-8"))
        checks.extend([
            {
                "name": "rendered_collision_audit",
                "passed": collision.get("verdict") == "PASS"
                and collision.get("summary", {}).get("fail") == 0,
                "verdict": collision.get("verdict"),
                "summary": collision.get("summary"),
            },
            {
                "name": "pdf_text_glyph_floor",
                "passed": pdf_text.get("auditable") is True
                and pdf_text.get("minimum_found_pt", 0) >= 5.0
                and pdf_text.get("below_minimum_count") == 0,
                "minimum_found_pt": pdf_text.get("minimum_found_pt"),
                "text_run_count": pdf_text.get("text_run_count"),
            },
            {
                "name": "single_composite_panel_alignment_gate",
                "passed": alignment.get("verdict") == "NOT APPLICABLE",
                "verdict": alignment.get("verdict"),
            },
        ])

    outputs = {}
    for path in [
        assembly_dir / "lan_host_review_filtered.fasta.gz",
        assembly_dir / "lan_removed_confirmed_contaminants.fasta.gz",
        plot_dir / "lan_review_filtered_blobplot.png",
        plot_dir / "lan_review_filtered_blobplot.svg",
        plot_dir / "lan_review_filtered_blobplot.pdf",
        audit_dir / "lan_final_decontam_review.tsv",
        audit_dir / "decontam_review_audit.json",
    ]:
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"Required output missing or empty: {path}")
        outputs[str(path.relative_to(root))] = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    failed = [check["name"] for check in checks if not check["passed"]]
    report = {
        "status": "passed" if not failed else "failed",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "expected": expected,
        "observed": observed,
        "checks": checks,
        "outputs": outputs,
        "failed_checks": failed,
    }
    (root / "final_qa.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if failed:
        raise RuntimeError("Validation failed: " + ", ".join(failed))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
