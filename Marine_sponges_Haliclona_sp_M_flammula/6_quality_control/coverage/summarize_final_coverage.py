#!/usr/bin/env python3
"""Audit and summarize matched ONT coverage for the two locked assemblies."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import statistics
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path("/work1/sponge/lht-njlan-final-coverage-20261009")
THRESHOLDS = (1, 5, 10, 20)
SCRIPT_VERSION = "1.0.1"


SAMPLES = OrderedDict(
    [
        (
            "LHT-04",
            {
                "fasta": Path(
                    "/work1/sponge/lht_reassembly_compare/15_read_validated_curated/"
                    "LHT_curated_primary.fasta"
                ),
                "reads": Path(
                    "/public/home/sponge/sci_data/lht/00_input/LHT04.ont.fq.gz"
                ),
                "out": ROOT / "LHT-04",
                "prefix": "LHT04.final844.ONT",
                "expected_contigs": 844,
                "expected_span": 115_734_075,
                "expected_fasta_sha256": (
                    "026d8b2e3e9888de0763d42ce226fcf7c63062a145e6ad8478e2e8c1725b2499"
                ),
                "expected_reads_sha256": (
                    "d75566c822e74e5807eb2b218de4bb8f585ef68a5600d6062bd7ecb4f9929a4b"
                ),
                "expected_bam_path": ROOT / "LHT-04" / "LHT04.final844.ONT.sorted.bam",
            },
        ),
        (
            "NJ-lan",
            {
                "fasta": Path(
                    "/work1/sponge/lan_tier2_decontam_review/assembly/"
                    "lan_host_review_filtered.fasta"
                ),
                "reads": Path(
                    "/public/home/sponge/sci_data/lan/00_input/Lanm.ont.fq.gz"
                ),
                "out": ROOT / "NJ-lan",
                "prefix": "NJ-lan.final7117.ONT",
                "expected_contigs": 7_117,
                "expected_span": 291_889_606,
                "expected_fasta_sha256": (
                    "cea9b34517164f02103111a3a1dfa7cfd65ec95a61c224da1901036a4cca120a"
                ),
                "expected_reads_sha256": (
                    "3453e91714f7cc1b9bf9d18fbaac6dc82eed51e4f300644d7464bf08a0d00587"
                ),
                "expected_bam_path": Path(
                    "/work1/sponge/lan_tier2_decontam_review/mapping/"
                    "lan.review_filtered.ONT.sorted.bam"
                ),
                "expected_bam_sha256": (
                    "d607aeac091f34ef4eca9e7cc9c091600de8bb857e99e1923713d2074b4cb949"
                ),
            },
        ),
    ]
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> OrderedDict[str, dict[str, int | float]]:
    records: OrderedDict[str, dict[str, int | float]] = OrderedDict()
    ident: str | None = None
    length = gc = acgt = n_bases = 0

    def save() -> None:
        nonlocal ident, length, gc, acgt, n_bases
        if ident is None:
            return
        if ident in records:
            raise ValueError(f"Duplicate FASTA identifier: {ident}")
        records[ident] = {
            "length_bp": length,
            "gc_bases": gc,
            "acgt_bases": acgt,
            "gc_fraction": gc / acgt if acgt else 0.0,
            "n_bases": n_bases,
        }

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                save()
                ident = line[1:].strip().split()[0]
                length = gc = acgt = n_bases = 0
                continue
            seq = line.strip().upper()
            length += len(seq)
            gc += seq.count("G") + seq.count("C")
            acgt += sum(seq.count(base) for base in "ACGT")
            n_bases += seq.count("N")
    save()
    if not records:
        raise ValueError(f"No FASTA records: {path}")
    return records


def read_seqkit(path: Path) -> dict[str, int | float | str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one seqkit row in {path}, got {len(rows)}")
    row = rows[0]
    if "N50" not in row or not row["N50"]:
        raise ValueError(
            f"Missing N50 in {path}; generate this table with seqkit stats -a -T"
        )
    return {
        "file": row["file"],
        "num_seqs": int(row["num_seqs"]),
        "sum_len": int(row["sum_len"]),
        "min_len": int(row["min_len"]),
        "avg_len": float(row["avg_len"]),
        "max_len": int(row["max_len"]),
        "n50": int(row["N50"]),
    }


def plus_total(match: re.Match[str]) -> int:
    return int(match.group(1)) + int(match.group(2))


def read_flagstat(path: Path) -> dict[str, int]:
    patterns = OrderedDict(
        [
            ("alignment_records", re.compile(r"^(\d+) \+ (\d+) in total ")),
            ("primary_records", re.compile(r"^(\d+) \+ (\d+) primary$")),
            ("secondary_records", re.compile(r"^(\d+) \+ (\d+) secondary$")),
            ("supplementary_records", re.compile(r"^(\d+) \+ (\d+) supplementary$")),
            ("mapped_alignment_records", re.compile(r"^(\d+) \+ (\d+) mapped \(")),
            ("primary_mapped_reads", re.compile(r"^(\d+) \+ (\d+) primary mapped \(")),
        ]
    )
    values: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        for key, pattern in patterns.items():
            match = pattern.match(line)
            if match:
                values[key] = plus_total(match)
    missing = set(patterns) - set(values)
    if missing:
        raise ValueError(f"Missing flagstat fields in {path}: {sorted(missing)}")
    return values


def read_mapping_counts(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["metric", "count"]:
            raise ValueError(f"Unexpected mapping-count header in {path}: {reader.fieldnames}")
        for row in reader:
            result[row["metric"]] = int(row["count"])
    required = {
        "primary_records_F2816",
        "primary_mapped_F2820",
        "primary_mapped_q20_F2820",
    }
    if not required.issubset(result):
        raise ValueError(f"Missing mapping counts in {path}: {sorted(required - set(result))}")
    return result


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_regions(path: Path) -> dict[str, dict[str, int | float]]:
    rows: dict[str, dict[str, int | float]] = {}
    with open_text(path) as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected >=4 columns")
            ident = fields[0]
            start, end, mean_depth = int(fields[1]), int(fields[2]), float(fields[3])
            if ident in rows:
                raise ValueError(f"Duplicate region: {ident}")
            if start != 0 or end <= 0 or mean_depth < 0 or not math.isfinite(mean_depth):
                raise ValueError(f"Invalid region in {path}:{line_number}")
            rows[ident] = {"length_bp": end, "mean_depth": mean_depth}
    return rows


def read_thresholds(path: Path) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"#chrom", "start", "end", *(f"{x}X" for x in THRESHOLDS)}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"Unexpected threshold header in {path}: {reader.fieldnames}")
        for line_number, row in enumerate(reader, 2):
            ident = row["#chrom"]
            if ident in result:
                raise ValueError(f"Duplicate threshold region: {ident}")
            start, end = int(row["start"]), int(row["end"])
            covered = {f"ge_{x}x_bp": int(row[f"{x}X"]) for x in THRESHOLDS}
            if start != 0 or end <= 0:
                raise ValueError(f"Invalid threshold interval in {path}:{line_number}")
            previous = end
            for threshold in THRESHOLDS:
                value = covered[f"ge_{threshold}x_bp"]
                if value < 0 or value > previous:
                    raise ValueError(
                        f"Non-monotonic threshold coverage for {ident} at {threshold}X"
                    )
                previous = value
            result[ident] = {"length_bp": end, **covered}
    return result


def read_mosdepth_total(path: Path) -> dict[str, int | float]:
    total: dict[str, int | float] | None = None
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "length", "bases", "mean", "min", "max"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected mosdepth summary header in {path}: {reader.fieldnames}")
        for row in reader:
            if row["chrom"] == "total":
                total = {
                    "reported_nonzero_reference_length": int(row["length"]),
                    "depth_bases": int(row["bases"]),
                    "reported_mean": float(row["mean"]),
                    "min_depth": int(row["min"]),
                    "max_depth": int(row["max"]),
                }
                break
    if total is None:
        raise ValueError(f"Missing total row in mosdepth summary: {path}")
    return total


def weighted_median(pairs: Iterable[tuple[float, int]]) -> float:
    ordered = sorted(pairs, key=lambda pair: pair[0])
    total_weight = sum(weight for _, weight in ordered)
    cumulative = 0
    for value, weight in ordered:
        cumulative += weight
        if cumulative * 2 >= total_weight:
            return value
    raise ValueError("Cannot calculate weighted median")


def n50(lengths: Iterable[int]) -> int:
    values = sorted(lengths, reverse=True)
    target = sum(values) / 2
    cumulative = 0
    for value in values:
        cumulative += value
        if cumulative >= target:
            return value
    raise ValueError("Cannot calculate N50")


def fmt_pct(value: float) -> str:
    return f"{value:.4f}%"


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_versions(out: Path) -> dict[str, str]:
    versions = {
        "minimap2": (out / "minimap2.version.txt").read_text(encoding="utf-8").splitlines()[0].strip(),
        "samtools": (out / "samtools.version.txt").read_text(encoding="utf-8").splitlines()[0].strip(),
        "mosdepth": (out / "mosdepth.version.txt").read_text(encoding="utf-8").splitlines()[0].strip(),
        "seqkit": (out / "seqkit.version.txt").read_text(encoding="utf-8").splitlines()[0].strip(),
        "python": (out / "python.version.txt").read_text(encoding="utf-8").splitlines()[0].strip(),
    }
    expected = {
        "minimap2": "2.31-r1302",
        "samtools": "samtools 1.23.1",
        "mosdepth": "mosdepth 0.3.14",
        "seqkit": "seqkit v2.13.0",
    }
    for tool, value in expected.items():
        if versions[tool] != value:
            raise ValueError(f"Unexpected {tool} version: {versions[tool]} != {value}")
    return versions


def main() -> int:
    summary_dir = Path(os.environ.get("COVERAGE_SUMMARY_DIR", str(ROOT / "summary")))
    summary_dir.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []
    warnings: list[str] = []
    all_metrics: list[dict[str, Any]] = []
    long_metrics: list[dict[str, Any]] = []
    sample_audits: dict[str, Any] = {}

    for sample, cfg in SAMPLES.items():
        out = Path(cfg["out"])
        prefix = str(cfg["prefix"])
        if not (out / "COMPLETED.txt").is_file():
            raise ValueError(f"Missing completion marker for {sample}")

        fasta = read_fasta(Path(cfg["fasta"]))
        fasta_sha = sha256(Path(cfg["fasta"]))
        contig_count = len(fasta)
        span = sum(int(row["length_bp"]) for row in fasta.values())
        if contig_count != cfg["expected_contigs"]:
            errors.append(f"{sample}: contig count {contig_count} != expected {cfg['expected_contigs']}")
        if span != cfg["expected_span"]:
            errors.append(f"{sample}: span {span} != expected {cfg['expected_span']}")
        if fasta_sha != cfg["expected_fasta_sha256"]:
            errors.append(f"{sample}: FASTA SHA-256 mismatch")

        seqkit = read_seqkit(out / "reads.seqkit.all.tsv")
        flagstat = read_flagstat(out / f"{prefix}.flagstat.txt")
        mapping_counts = read_mapping_counts(out / "mapping_counts.tsv")
        versions = read_versions(out)
        if Path(str(seqkit["file"])) != Path(cfg["reads"]):
            errors.append(f"{sample}: seqkit input path does not equal the locked reads path")
        if seqkit["num_seqs"] != flagstat["primary_records"]:
            errors.append(
                f"{sample}: FASTQ reads {seqkit['num_seqs']} != flagstat primary records "
                f"{flagstat['primary_records']}"
            )
        if flagstat["alignment_records"] != (
            flagstat["primary_records"]
            + flagstat["secondary_records"]
            + flagstat["supplementary_records"]
        ):
            errors.append(f"{sample}: flagstat record categories do not sum to total")
        if not (0 < flagstat["primary_mapped_reads"] <= flagstat["primary_records"]):
            errors.append(f"{sample}: invalid primary mapped/total counts")
        if not (0 <= flagstat["mapped_alignment_records"] <= flagstat["alignment_records"]):
            errors.append(f"{sample}: invalid mapped alignment-record count")
        if mapping_counts["primary_records_F2816"] != flagstat["primary_records"]:
            errors.append(f"{sample}: samtools view primary denominator != flagstat primary")
        if mapping_counts["primary_mapped_F2820"] != flagstat["primary_mapped_reads"]:
            errors.append(f"{sample}: samtools view primary mapped != flagstat primary mapped")
        if not (
            0
            <= mapping_counts["primary_mapped_q20_F2820"]
            <= mapping_counts["primary_mapped_F2820"]
        ):
            errors.append(f"{sample}: invalid MAPQ>=20 primary mapped count")

        reads_sha_line = (out / "reads.sha256").read_text(encoding="utf-8").strip().split()
        if (
            len(reads_sha_line) != 2
            or re.fullmatch(r"[0-9a-f]{64}", reads_sha_line[0]) is None
            or Path(reads_sha_line[1]) != Path(cfg["reads"])
        ):
            errors.append(f"{sample}: invalid reads checksum record")
            reads_sha = ""
        else:
            reads_sha = reads_sha_line[0]
            if reads_sha != cfg["expected_reads_sha256"]:
                errors.append(f"{sample}: reads SHA-256 mismatch")

        bam_sha_line = (out / "bam.sha256").read_text(encoding="utf-8").splitlines()[0].split()
        if (
            len(bam_sha_line) != 2
            or re.fullmatch(r"[0-9a-f]{64}", bam_sha_line[0]) is None
            or Path(bam_sha_line[1]) != Path(cfg["expected_bam_path"])
        ):
            errors.append(f"{sample}: invalid BAM checksum record")
            bam_sha = ""
        else:
            bam_sha = bam_sha_line[0]
            if "expected_bam_sha256" in cfg and bam_sha != cfg["expected_bam_sha256"]:
                errors.append(f"{sample}: BAM SHA-256 mismatch")

        paths = {
            "default_regions": out / f"{prefix}.default.regions.bed.gz",
            "default_thresholds": out / f"{prefix}.default.thresholds.bed.gz",
            "default_summary": out / f"{prefix}.default.mosdepth.summary.txt",
            "primary_regions": out / f"{prefix}.primary.regions.bed.gz",
            "primary_thresholds": out / f"{prefix}.primary.thresholds.bed.gz",
            "primary_summary": out / f"{prefix}.primary.mosdepth.summary.txt",
            "primary_q20_regions": out / f"{prefix}.primary_q20.regions.bed.gz",
            "primary_q20_thresholds": out / f"{prefix}.primary_q20.thresholds.bed.gz",
            "primary_q20_summary": out / f"{prefix}.primary_q20.mosdepth.summary.txt",
        }
        for path in paths.values():
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing coverage result: {path}")

        default_regions = read_regions(paths["default_regions"])
        primary_regions = read_regions(paths["primary_regions"])
        primary_q20_regions = read_regions(paths["primary_q20_regions"])
        default_thresholds = read_thresholds(paths["default_thresholds"])
        primary_thresholds = read_thresholds(paths["primary_thresholds"])
        primary_q20_thresholds = read_thresholds(paths["primary_q20_thresholds"])
        depth_totals = {
            "default": read_mosdepth_total(paths["default_summary"]),
            "primary": read_mosdepth_total(paths["primary_summary"]),
            "primary_q20": read_mosdepth_total(paths["primary_q20_summary"]),
        }
        expected_ids = set(fasta)
        for label, records in (
            ("default regions", default_regions),
            ("primary regions", primary_regions),
            ("primary Q20 regions", primary_q20_regions),
            ("default thresholds", default_thresholds),
            ("primary thresholds", primary_thresholds),
            ("primary Q20 thresholds", primary_q20_thresholds),
        ):
            if set(records) != expected_ids:
                raise ValueError(
                    f"{sample}: {label} ID mismatch: FASTA-only="
                    f"{len(expected_ids - set(records))}, result-only="
                    f"{len(set(records) - expected_ids)}"
                )
            length_mismatches = [
                ident
                for ident in expected_ids & set(records)
                if int(records[ident]["length_bp"]) != int(fasta[ident]["length_bp"])
            ]
            if length_mismatches:
                raise ValueError(f"{sample}: {label} has {len(length_mismatches)} length mismatches")

        per_contig_rows: list[dict[str, Any]] = []
        for ident, sequence in fasta.items():
            length = int(sequence["length_bp"])
            row: dict[str, Any] = {
                "sample_id": sample,
                "contig_id": ident,
                "length_bp": length,
                "gc_fraction": float(sequence["gc_fraction"]),
                "n_bases": int(sequence["n_bases"]),
                "default_mean_depth_x": float(default_regions[ident]["mean_depth"]),
                "primary_mean_depth_x": float(primary_regions[ident]["mean_depth"]),
                "primary_q20_mean_depth_x": float(primary_q20_regions[ident]["mean_depth"]),
            }
            for label, thresholds in (
                ("default", default_thresholds),
                ("primary", primary_thresholds),
                ("primary_q20", primary_q20_thresholds),
            ):
                for threshold in THRESHOLDS:
                    covered = thresholds[ident][f"ge_{threshold}x_bp"]
                    row[f"{label}_covered_bp_ge_{threshold}x"] = covered
                    row[f"{label}_breadth_pct_ge_{threshold}x"] = covered / length * 100
            per_contig_rows.append(row)

        per_contig_fields = list(per_contig_rows[0])
        write_tsv(summary_dir / f"{sample}.per_contig_coverage.tsv", per_contig_fields, per_contig_rows)

        primary_depth = [float(row["primary_mean_depth_x"]) for row in per_contig_rows]
        primary_q20_depth = [float(row["primary_q20_mean_depth_x"]) for row in per_contig_rows]
        default_depth = [float(row["default_mean_depth_x"]) for row in per_contig_rows]
        metric: dict[str, Any] = {
            "sample_id": sample,
            "assembly_contigs": contig_count,
            "assembly_span_bp": span,
            "assembly_n50_bp": n50(int(row["length_bp"]) for row in fasta.values()),
            "assembly_gc_percent": (
                sum(int(row["gc_bases"]) for row in fasta.values())
                / sum(int(row["acgt_bases"]) for row in fasta.values())
                * 100
            ),
            "raw_ont_reads": seqkit["num_seqs"],
            "raw_ont_bases": seqkit["sum_len"],
            "raw_ont_yield_gb": float(seqkit["sum_len"]) / 1_000_000_000,
            "raw_ont_read_n50_bp": seqkit["n50"],
            "primary_mapped_reads": flagstat["primary_mapped_reads"],
            "primary_unmapped_reads": flagstat["primary_records"] - flagstat["primary_mapped_reads"],
            "primary_mapping_percent": flagstat["primary_mapped_reads"] / flagstat["primary_records"] * 100,
            "primary_mapped_q20_reads": mapping_counts["primary_mapped_q20_F2820"],
            "primary_mapping_q20_percent": mapping_counts["primary_mapped_q20_F2820"] / flagstat["primary_records"] * 100,
            "alignment_records": flagstat["alignment_records"],
            "mapped_alignment_records": flagstat["mapped_alignment_records"],
            "mapped_alignment_record_percent": flagstat["mapped_alignment_records"] / flagstat["alignment_records"] * 100,
            "secondary_alignment_records": flagstat["secondary_records"],
            "supplementary_alignment_records": flagstat["supplementary_records"],
            "default_depth_bases": depth_totals["default"]["depth_bases"],
            "default_length_weighted_mean_depth_x": int(depth_totals["default"]["depth_bases"]) / span,
            "primary_depth_bases": depth_totals["primary"]["depth_bases"],
            "primary_length_weighted_mean_depth_x": int(depth_totals["primary"]["depth_bases"]) / span,
            "primary_q20_depth_bases": depth_totals["primary_q20"]["depth_bases"],
            "primary_q20_length_weighted_mean_depth_x": int(depth_totals["primary_q20"]["depth_bases"]) / span,
            "primary_contig_mean_depth_median_x": statistics.median(primary_depth),
            "primary_length_weighted_median_contig_mean_depth_x": weighted_median(
                (float(row["primary_mean_depth_x"]), int(row["length_bp"]))
                for row in per_contig_rows
            ),
            "primary_q20_contig_mean_depth_median_x": statistics.median(primary_q20_depth),
            "primary_q20_length_weighted_median_contig_mean_depth_x": weighted_median(
                (float(row["primary_q20_mean_depth_x"]), int(row["length_bp"]))
                for row in per_contig_rows
            ),
            "fasta_sha256": fasta_sha,
            "reads_sha256": reads_sha,
            "bam_sha256": bam_sha,
        }
        for label in ("default", "primary", "primary_q20"):
            for threshold in THRESHOLDS:
                key = f"{label}_covered_bp_ge_{threshold}x"
                covered = sum(int(row[key]) for row in per_contig_rows)
                metric[key] = covered
                metric[f"{label}_breadth_percent_ge_{threshold}x"] = covered / span * 100
                metric[f"{label}_zero_breadth_contigs_ge_{threshold}x"] = sum(
                    int(row[key]) == 0 for row in per_contig_rows
                )
            metric[f"{label}_zero_coverage_bp"] = span - int(metric[f"{label}_covered_bp_ge_1x"])
            metric[f"{label}_zero_coverage_percent"] = metric[f"{label}_zero_coverage_bp"] / span * 100
            metric[f"{label}_fully_zero_contigs"] = int(metric[f"{label}_zero_breadth_contigs_ge_1x"])
            metric[f"{label}_fully_zero_contig_span_bp"] = sum(
                int(row["length_bp"])
                for row in per_contig_rows
                if int(row[f"{label}_covered_bp_ge_1x"]) == 0
            )
            if not (
                0
                <= int(metric[f"{label}_fully_zero_contig_span_bp"])
                <= int(metric[f"{label}_zero_coverage_bp"])
                <= span
            ):
                errors.append(f"{sample}: invalid zero-coverage accounting for {label}")

        if not (
            int(metric["primary_q20_depth_bases"])
            <= int(metric["primary_depth_bases"])
            <= int(metric["default_depth_bases"])
        ):
            errors.append(f"{sample}: depth totals violate primary_q20 <= primary <= default")
        for threshold in THRESHOLDS:
            if not (
                int(metric[f"primary_q20_covered_bp_ge_{threshold}x"])
                <= int(metric[f"primary_covered_bp_ge_{threshold}x"])
                <= int(metric[f"default_covered_bp_ge_{threshold}x"])
            ):
                errors.append(
                    f"{sample}: {threshold}X breadth violates primary_q20 <= primary <= default"
                )
        all_metrics.append(metric)

        article_spec = [
            ("assembly_contigs", "Assembly contigs", "contigs", "Locked final FASTA"),
            ("assembly_span_bp", "Assembly span", "bp", "Locked final FASTA"),
            ("raw_ont_reads", "Raw ONT reads", "reads", "seqkit stats"),
            ("raw_ont_bases", "Raw ONT yield", "bp", "seqkit stats"),
            ("primary_mapped_reads", "Primary mapped ONT reads", "reads", "samtools flagstat"),
            ("primary_mapping_percent", "Primary ONT mapping rate", "%", "primary mapped / raw reads"),
            ("primary_mapped_q20_reads", "Primary mapped ONT reads with MAPQ >=20", "reads", "samtools view -F 2820 -q 20"),
            ("primary_mapping_q20_percent", "Primary ONT mapping rate with MAPQ >=20", "%", "primary MAPQ>=20 mapped / raw reads"),
            ("default_length_weighted_mean_depth_x", "Compatibility ONT depth", "x", "mosdepth -F 1796 -Q 0; includes supplementary alignments"),
            ("primary_q20_length_weighted_mean_depth_x", "High-confidence primary-only ONT depth", "x", "mosdepth -F 3844 -Q 20"),
            ("primary_q20_covered_bp_ge_1x", "High-confidence assembly covered at >=1x", "bp", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_breadth_percent_ge_1x", "High-confidence assembly breadth at >=1x", "%", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_covered_bp_ge_5x", "High-confidence assembly covered at >=5x", "bp", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_breadth_percent_ge_5x", "High-confidence assembly breadth at >=5x", "%", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_covered_bp_ge_10x", "High-confidence assembly covered at >=10x", "bp", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_breadth_percent_ge_10x", "High-confidence assembly breadth at >=10x", "%", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_covered_bp_ge_20x", "High-confidence assembly covered at >=20x", "bp", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_breadth_percent_ge_20x", "High-confidence assembly breadth at >=20x", "%", "mosdepth -F 3844 -Q 20 thresholds"),
            ("primary_q20_zero_coverage_bp", "Bases with zero high-confidence primary ONT coverage", "bp", "assembly span - high-confidence bases at >=1x"),
            ("primary_q20_zero_coverage_percent", "Assembly with zero high-confidence primary ONT coverage", "%", "assembly span - high-confidence bases at >=1x"),
            ("primary_q20_fully_zero_contigs", "Fully zero-covered contigs", "contigs", "high-confidence bases at >=1x = 0"),
            ("primary_q20_fully_zero_contig_span_bp", "Span of fully zero-covered contigs", "bp", "high-confidence bases at >=1x = 0"),
        ]
        for key, label, unit, method in article_spec:
            long_metrics.append(
                {
                    "sample_id": sample,
                    "metric_key": key,
                    "metric": label,
                    "value": metric[key],
                    "unit": unit,
                    "method": method,
                }
            )

        sample_audits[sample] = {
            "fasta": str(cfg["fasta"]),
            "reads": str(cfg["reads"]),
            "fasta_sha256": fasta_sha,
            "reads_sha256": reads_sha,
            "contigs": contig_count,
            "span_bp": span,
            "fastq_records_equal_flagstat_primary_records": seqkit["num_seqs"] == flagstat["primary_records"],
            "coverage_id_sets_equal_fasta": all(
                set(records) == expected_ids
                for records in (
                    default_regions,
                    primary_regions,
                    primary_q20_regions,
                    default_thresholds,
                    primary_thresholds,
                    primary_q20_thresholds,
                )
            ),
            "bam_sha256": bam_sha,
            "versions": versions,
            "coverage_definitions": {
                "default_compatibility": "mosdepth -F 1796 -Q 0; supplementary alignments included",
                "primary_q0": "mosdepth -F 3844 -Q 0; supplementary alignments excluded",
                "high_confidence_primary": "mosdepth -F 3844 -Q 20; supplementary alignments excluded and MAPQ >=20 required",
            },
        }

    if errors:
        status = "FAIL"
    else:
        status = "PASS"
    warnings.append(
        "The compatibility depth uses mosdepth's traditional exclusion mask 1796 and includes supplementary alignments; it is retained for continuity with the blob plots."
    )
    warnings.append(
        "Whole-assembly mean depth was calculated from integer mosdepth depth-bases divided by the complete assembly span; per-contig means remain rounded to two decimal places."
    )

    metrics_fields = list(all_metrics[0])
    write_tsv(summary_dir / "final_coverage_metrics.tsv", metrics_fields, all_metrics)
    write_tsv(
        summary_dir / "article_coverage_metrics_long.tsv",
        ["sample_id", "metric_key", "metric", "value", "unit", "method"],
        long_metrics,
    )

    audit = {
        "script_version": SCRIPT_VERSION,
        "status": status,
        "errors": errors,
        "warnings": warnings,
        "sample_audits": sample_audits,
    }
    (summary_dir / "coverage_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    if errors:
        raise ValueError("; ".join(errors))

    by_sample = {row["sample_id"]: row for row in all_metrics}
    lht = by_sample["LHT-04"]
    lan = by_sample["NJ-lan"]
    readme = f"""# LHT-04 与 NJ-lan 最终组装 ONT 覆盖结果

本目录给出锁定最终组装上的一致口径 ONT 重比对和覆盖统计。LHT-04 使用最终 844-contig FASTA 重新比对完整原始 ONT reads；NJ-lan 复核已完成的 7,117-contig 重比对 BAM，并补算相同的 primary-only 覆盖指标。

## 文章建议采用的结果

| 指标 | LHT-04 | NJ-lan |
|---|---:|---:|
| 最终组装 | {lht['assembly_contigs']:,} contigs / {lht['assembly_span_bp']:,} bp | {lan['assembly_contigs']:,} contigs / {lan['assembly_span_bp']:,} bp |
| 原始 ONT reads | {lht['raw_ont_reads']:,} / {lht['raw_ont_bases']:,} bp | {lan['raw_ont_reads']:,} / {lan['raw_ont_bases']:,} bp |
| Primary mapped reads | {lht['primary_mapped_reads']:,} ({fmt_pct(lht['primary_mapping_percent'])}) | {lan['primary_mapped_reads']:,} ({fmt_pct(lan['primary_mapping_percent'])}) |
| Primary mapped reads，MAPQ >=20 | {lht['primary_mapped_q20_reads']:,} ({fmt_pct(lht['primary_mapping_q20_percent'])}) | {lan['primary_mapped_q20_reads']:,} ({fmt_pct(lan['primary_mapping_q20_percent'])}) |
| 兼容旧 blob 图的平均深度（F1796/Q0） | {lht['default_length_weighted_mean_depth_x']:.2f}x | {lan['default_length_weighted_mean_depth_x']:.2f}x |
| 高置信 primary-only 平均深度（F3844/Q20） | {lht['primary_q20_length_weighted_mean_depth_x']:.2f}x | {lan['primary_q20_length_weighted_mean_depth_x']:.2f}x |
| 高置信覆盖广度 >=1x | {lht['primary_q20_covered_bp_ge_1x']:,} bp ({fmt_pct(lht['primary_q20_breadth_percent_ge_1x'])}) | {lan['primary_q20_covered_bp_ge_1x']:,} bp ({fmt_pct(lan['primary_q20_breadth_percent_ge_1x'])}) |
| 高置信覆盖广度 >=5x | {lht['primary_q20_covered_bp_ge_5x']:,} bp ({fmt_pct(lht['primary_q20_breadth_percent_ge_5x'])}) | {lan['primary_q20_covered_bp_ge_5x']:,} bp ({fmt_pct(lan['primary_q20_breadth_percent_ge_5x'])}) |
| 高置信覆盖广度 >=10x | {lht['primary_q20_covered_bp_ge_10x']:,} bp ({fmt_pct(lht['primary_q20_breadth_percent_ge_10x'])}) | {lan['primary_q20_covered_bp_ge_10x']:,} bp ({fmt_pct(lan['primary_q20_breadth_percent_ge_10x'])}) |
| 高置信覆盖广度 >=20x | {lht['primary_q20_covered_bp_ge_20x']:,} bp ({fmt_pct(lht['primary_q20_breadth_percent_ge_20x'])}) | {lan['primary_q20_covered_bp_ge_20x']:,} bp ({fmt_pct(lan['primary_q20_breadth_percent_ge_20x'])}) |
| 高置信零覆盖碱基 | {lht['primary_q20_zero_coverage_bp']:,} bp ({fmt_pct(lht['primary_q20_zero_coverage_percent'])}) | {lan['primary_q20_zero_coverage_bp']:,} bp ({fmt_pct(lan['primary_q20_zero_coverage_percent'])}) |
| 完全零覆盖 contig | {lht['primary_q20_fully_zero_contigs']:,} / {lht['primary_q20_fully_zero_contig_span_bp']:,} bp | {lan['primary_q20_fully_zero_contigs']:,} / {lan['primary_q20_fully_zero_contig_span_bp']:,} bp |

## 统计定义

- Primary mapping rate = `samtools flagstat` 的 primary mapped records / primary records；分母经核对等于 FASTQ read 数。
- 可与既有 NJ-lan blob 图直接比较的 canonical 深度使用 mosdepth exclusion mask 1796、MAPQ 0；该口径包含 supplementary alignments。
- 高置信技术验证使用 exclusion mask 3844、MAPQ >=20，排除 unmapped、secondary、supplementary、QC-failed 和 duplicate alignments。
- 全组装平均深度按 contig 长度加权；覆盖广度为达到 1x、5x、10x 或 20x 的碱基数除以完整组装长度。
- primary-only Q0 指标也保留在完整 TSV 中，用来分辨 supplementary 过滤与 MAPQ 过滤各自的影响。
- NJ-lan 早期 7,366-contig blob 结果（51.45x）不是本稿锁定的 7,117-contig 最终组装结果，不能用于文章。

## QA

- 最终 FASTA 的 contig 数、总长和 SHA-256 与锁定记录完全一致。
- FASTQ read 数与 BAM primary record 数完全一致。
- BAM reference ID/长度、mosdepth region ID/长度和 FASTA 逐条一致。
- 1x/5x/10x/20x 覆盖碱基数逐条单调且不超过 contig 长度。
- 最终审计状态：PASS。
"""
    (summary_dir / "README_结果说明.md").write_text(readme, encoding="utf-8")

    methods = (
        "Oxford Nanopore reads were aligned to each locked final assembly with "
        "minimap2 v2.31-r1302 using the map-ont preset (-ax map-ont -I 20G). "
        "Alignments were coordinate-sorted and indexed with SAMtools v1.23.1. "
        "Primary read mapping rates were calculated from samtools flagstat as the "
        "number of primary mapped records divided by the number of primary records, "
        "which was verified against the number of input FASTQ records. Coverage was "
        "calculated with mosdepth v0.3.14 over one full-length interval per contig. "
        "For direct continuity with the GC-coverage analysis, compatibility depth used "
        "mosdepth flag mask 1796 and MAPQ 0. A high-confidence primary-only analysis "
        "additionally excluded supplementary alignments and required MAPQ >=20 (flag "
        "mask 3844). Genome-wide mean depth was calculated from the integer sum of "
        "depth-bases divided by the complete assembly span, and high-confidence "
        "coverage breadth was reported at >=1x, >=5x, >=10x and >=20x.\n"
    )
    results = (
        f"For LHT-04, {lht['primary_mapped_reads']:,} of {lht['raw_ont_reads']:,} "
        f"ONT reads mapped as primary alignments ({lht['primary_mapping_percent']:.2f}%), "
        f"yielding a compatibility depth of {lht['default_length_weighted_mean_depth_x']:.2f}x "
        f"and a high-confidence primary-only depth of {lht['primary_q20_length_weighted_mean_depth_x']:.2f}x; "
        f"{lht['primary_q20_breadth_percent_ge_1x']:.2f}% of the assembly had high-confidence coverage at >=1x, "
        f"and {lht['primary_q20_fully_zero_contigs']:,} contigs were completely uncovered. "
        f"For NJ-lan, {lan['primary_mapped_reads']:,} of {lan['raw_ont_reads']:,} "
        f"ONT reads mapped as primary alignments ({lan['primary_mapping_percent']:.2f}%), "
        f"with a compatibility depth of {lan['default_length_weighted_mean_depth_x']:.2f}x "
        f"and a high-confidence primary-only depth of {lan['primary_q20_length_weighted_mean_depth_x']:.2f}x; "
        f"{lan['primary_q20_breadth_percent_ge_1x']:.2f}% of the assembly had high-confidence coverage at >=1x, "
        f"and {lan['primary_q20_fully_zero_contigs']:,} contigs were completely uncovered.\n"
    )
    (summary_dir / "article_methods_and_results.txt").write_text(
        "METHODS\n" + methods + "\nRESULTS\n" + results, encoding="utf-8"
    )
    print(json.dumps({"status": status, "metrics": all_metrics}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
