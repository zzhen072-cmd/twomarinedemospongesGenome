#!/usr/bin/env python3
"""Read-only audit of the final LHT-04 and NJ-lan annotation workflow.

The workflow inputs and result layout are intentionally fixed.  This script reads
those inputs/results and writes only the documented summary products under
``ROOT/08_summary``, ``ROOT/{sample}/08_summary``, and the four standardized
top-level summary/audit paths.  It does not repair or rename scientific outputs.
A failed check is retained in the JSON/TSV audit and causes a non-zero exit status.

Python >= 3.9; standard library only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import re
import socket
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, TextIO, Tuple
from urllib.parse import unquote


AUDIT_VERSION = "1.0.2"
EXPECTED_BUSCO_VERSION = "6.0.0"
EXPECTED_BUSCO_LINEAGE = "metazoa_odb10"
EXPECTED_BUSCO_GROUPS = 954
STAGES = (
    "input",
    "repeatmasker",
    "genome_busco",
    "braker",
    "protein_busco",
    "eggnog",
    "swissprot",
    "interpro",
)

SAMPLES = {
    "lht04": {
        "label": "LHT-04",
        "source": "/work1/sponge/lht_reassembly_compare/15_read_validated_curated/LHT_curated_primary.fasta",
        "records": 844,
        "bp": 115_734_075,
        "sha256": "026d8b2e3e9888de0763d42ce226fcf7c63062a145e6ad8478e2e8c1725b2499",
        "short_prefix": "lhtf_ctg",
        "braker_min_contig": 10_000,
    },
    "njlan": {
        "label": "NJ-lan",
        "source": "/work1/sponge/lan_tier2_decontam_review/assembly/lan_host_review_filtered.fasta",
        "records": 7_117,
        "bp": 291_889_606,
        "sha256": "cea9b34517164f02103111a3a1dfa7cfd65ec95a61c224da1901036a4cca120a",
        "short_prefix": "njlf_ctg",
        "braker_min_contig": 2_000,
    },
}

DNA_ALPHABET = set("ACGTNRYKMSWBDHV")
PROTEIN_ALPHABET = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")
MISSING_VALUES = {"", "-", ".", "NA", "N/A", "None", "null"}
IPR_RE = re.compile(r"^IPR\d{6}$")


@dataclass
class Issue:
    severity: str
    sample: str
    stage: str
    code: str
    message: str
    path: str = ""
    details: object = None


@dataclass
class FileRecord:
    sample: str
    stage: str
    role: str
    path: str
    bytes_before: int
    mtime_ns_before: int
    sha256: str = ""
    stable_during_audit: Optional[bool] = None


@dataclass
class FastaRecord:
    identifier: str
    header: str
    sequence: str


class Audit:
    """Collect checks without silently discarding later-stage evidence."""

    def __init__(self, root: Path, max_issues: int = 2000) -> None:
        self.root = root
        self.max_issues = max_issues
        self.issues: List[Issue] = []
        self.issue_totals: Counter = Counter()
        self.issue_totals_by_stage: Counter = Counter()
        self.suppressed_issue_totals: Counter = Counter()
        self.metrics: Dict[str, Dict[str, object]] = {sample: {} for sample in SAMPLES}
        self.busco_rows: List[dict] = []
        self.functional_rows: List[dict] = []
        self.files: Dict[str, FileRecord] = {}
        self.attempted_stages: Set[Tuple[str, str]] = set()

    def issue(
        self,
        severity: str,
        sample: str,
        stage: str,
        code: str,
        message: str,
        path: Optional[Path] = None,
        details: object = None,
    ) -> None:
        severity = severity.upper()
        self.issue_totals[severity] += 1
        self.issue_totals_by_stage[(sample, stage, severity)] += 1
        if len(self.issues) < self.max_issues:
            self.issues.append(
                Issue(
                    severity=severity,
                    sample=sample,
                    stage=stage,
                    code=code,
                    message=message,
                    path=str(path) if path else "",
                    details=details,
                )
            )
        else:
            self.suppressed_issue_totals[severity] += 1

    def error(self, *args, **kwargs) -> None:
        self.issue("ERROR", *args, **kwargs)

    def warning(self, *args, **kwargs) -> None:
        self.issue("WARNING", *args, **kwargs)

    def info(self, *args, **kwargs) -> None:
        self.issue("INFO", *args, **kwargs)

    def expect(
        self,
        condition: bool,
        sample: str,
        stage: str,
        code: str,
        message: str,
        path: Optional[Path] = None,
        details: object = None,
    ) -> bool:
        if not condition:
            self.error(sample, stage, code, message, path, details)
        return condition

    def add_file(self, path: Path, sample: str, stage: str, role: str) -> None:
        key = str(path)
        try:
            stat = path.stat()
        except OSError:
            return
        if key not in self.files:
            self.files[key] = FileRecord(
                sample=sample,
                stage=stage,
                role=role,
                path=key,
                bytes_before=stat.st_size,
                mtime_ns_before=stat.st_mtime_ns,
            )

    def required_file(
        self,
        path: Path,
        sample: str,
        stage: str,
        role: str,
        nonempty: bool = True,
    ) -> bool:
        if not path.is_file():
            self.error(sample, stage, "MISSING_FILE", f"Missing required {role}", path)
            return False
        if nonempty and path.stat().st_size == 0:
            self.error(sample, stage, "EMPTY_FILE", f"Required {role} is empty", path)
            return False
        self.add_file(path, sample, stage, role)
        return True

    def run_stage(self, sample: str, stage: str, function) -> None:
        self.attempted_stages.add((sample, stage))
        try:
            function()
        except Exception as exc:  # preserve a machine-readable report on malformed input
            self.error(
                sample,
                stage,
                "UNHANDLED_STAGE_EXCEPTION",
                f"{type(exc).__name__}: {exc}",
                details={"traceback": traceback.format_exc(limit=12)},
            )

    def finalize_file_hashes(self) -> None:
        """Hash each audited primary file and detect concurrent modification."""
        for record in self.files.values():
            path = Path(record.path)
            try:
                record.sha256 = sha256_file(path)
                after = path.stat()
                record.stable_during_audit = (
                    after.st_size == record.bytes_before
                    and after.st_mtime_ns == record.mtime_ns_before
                )
                if not record.stable_during_audit:
                    self.error(
                        record.sample,
                        record.stage,
                        "FILE_CHANGED_DURING_AUDIT",
                        "File size or modification time changed while it was audited",
                        path,
                        {
                            "bytes_before": record.bytes_before,
                            "bytes_after": after.st_size,
                            "mtime_ns_before": record.mtime_ns_before,
                            "mtime_ns_after": after.st_mtime_ns,
                        },
                    )
            except Exception as exc:
                record.stable_during_audit = False
                self.error(
                    record.sample,
                    record.stage,
                    "FILE_HASH_FAILED",
                    f"Could not hash audited file: {exc}",
                    path,
                )


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("rt", encoding="utf-8", errors="strict")


def open_uncompressed_binary(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_uncompressed(path: Path) -> str:
    digest = hashlib.sha256()
    with open_uncompressed_binary(path) as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_integrity(path: Path) -> None:
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(1024 * 1024), b""):
            pass


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
    identifier: Optional[str] = None
    header = ""
    sequence: List[str] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if identifier is not None:
                    yield FastaRecord(identifier, header, "".join(sequence))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                identifier = header.split()[0]
                sequence = []
            else:
                if identifier is None:
                    if line.strip():
                        raise ValueError(
                            f"{path}:{line_number}: sequence data before first header"
                        )
                    continue
                sequence.append("".join(line.split()))
    if identifier is not None:
        yield FastaRecord(identifier, header, "".join(sequence))


def read_fasta_map(
    path: Path,
    audit: Audit,
    sample: str,
    stage: str,
    alphabet: Set[str],
    sequence_kind: str,
) -> Tuple[Dict[str, str], List[str]]:
    records: Dict[str, str] = {}
    order: List[str] = []
    for record in iter_fasta(path):
        if record.identifier in records:
            audit.error(
                sample,
                stage,
                "DUPLICATE_FASTA_ID",
                f"Duplicate {sequence_kind} FASTA identifier {record.identifier}",
                path,
            )
            continue
        if not record.sequence:
            audit.error(
                sample,
                stage,
                "EMPTY_FASTA_SEQUENCE",
                f"Empty {sequence_kind} sequence {record.identifier}",
                path,
            )
        invalid = sorted(set(record.sequence.upper()) - alphabet)
        if invalid:
            audit.error(
                sample,
                stage,
                "INVALID_FASTA_SYMBOLS",
                f"Unsupported {sequence_kind} symbols in {record.identifier}: {''.join(invalid)}",
                path,
            )
        records[record.identifier] = record.sequence.upper()
        order.append(record.identifier)
    if not records:
        audit.error(
            sample,
            stage,
            "NO_FASTA_RECORDS",
            f"No {sequence_kind} records found",
            path,
        )
    return records, order


def parse_json_file(path: Path) -> dict:
    with path.open("rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def compare_sets(
    audit: Audit,
    sample: str,
    stage: str,
    code: str,
    left_name: str,
    left: Set[str],
    right_name: str,
    right: Set[str],
    path: Optional[Path] = None,
) -> bool:
    if left == right:
        return True
    audit.error(
        sample,
        stage,
        code,
        f"{left_name} and {right_name} differ",
        path,
        {
            f"only_{left_name}": sorted(left - right)[:20],
            f"only_{right_name}": sorted(right - left)[:20],
            f"{left_name}_count": len(left),
            f"{right_name}_count": len(right),
        },
    )
    return False


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("wt", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(str(temporary), str(path))


def tsv_text(rows: Sequence[dict], columns: Sequence[str]) -> str:
    from io import StringIO

    output = StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=list(columns),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        rendered = {}
        for key, value in row.items():
            if isinstance(value, set):
                value = sorted(value)
            rendered[key] = (
                json.dumps(value, ensure_ascii=False, sort_keys=True)
                if isinstance(value, (dict, list, tuple))
                else value
            )
        writer.writerow(
            rendered
        )
    return output.getvalue()


class SampleAuditor:
    def __init__(self, audit: Audit, sample: str, config: dict) -> None:
        self.audit = audit
        self.sample = sample
        self.config = config
        self.base = audit.root / sample
        self.metrics = audit.metrics[sample]
        self.short_lengths: Dict[str, int] = {}
        self.short_order: List[str] = []
        self.masked_lengths: Dict[str, int] = {}
        self.gff_gene_ids: Set[str] = set()
        self.gff_transcript_ids: Set[str] = set()
        self.tx_to_gene: Dict[str, str] = {}
        self.braker_proteins: Dict[str, str] = {}
        self.braker_protein_order: List[str] = []
        self.braker_cds: Dict[str, str] = {}
        self.representatives: Dict[str, str] = {}
        self.representative_order: List[str] = []

    def run(self) -> None:
        for stage, function in (
            ("input", self.audit_input),
            ("repeatmasker", self.audit_repeatmasker),
            ("genome_busco", self.audit_genome_busco),
            ("braker", self.audit_braker),
            ("protein_busco", self.audit_protein_and_busco),
            ("eggnog", self.audit_eggnog),
            ("swissprot", self.audit_swissprot),
            ("interpro", self.audit_interpro),
        ):
            self.audit.run_stage(self.sample, stage, function)

    def require(self, path: Path, stage: str, role: str, nonempty: bool = True) -> bool:
        return self.audit.required_file(
            path, self.sample, stage, role, nonempty=nonempty
        )

    def audit_input(self) -> None:
        stage = "input"
        root_audit_path = self.audit.root / "input_audit.json"
        marker = self.audit.root / ".inputs_prepared"
        manifest_path = self.base / "00_input" / "input_manifest.json"
        map_path = self.base / "00_input" / "contig_id_map.tsv"
        short_path = self.base / "00_input" / f"{self.sample}_final.shortids.fasta"
        required = [
            self.require(root_audit_path, stage, "combined input audit"),
            self.require(marker, stage, "input preparation marker"),
            self.require(manifest_path, stage, "sample input manifest"),
            self.require(map_path, stage, "contig ID map"),
            self.require(short_path, stage, "short-ID input FASTA"),
        ]
        if not all(required):
            return

        if marker.read_text(encoding="ascii").strip() != "PASS":
            self.audit.error(
                self.sample,
                stage,
                "INVALID_PREPARED_MARKER",
                "Input preparation marker is not PASS",
                marker,
            )

        combined = parse_json_file(root_audit_path)
        manifest = parse_json_file(manifest_path)
        self.audit.expect(
            combined.get("status") == "PASS",
            self.sample,
            stage,
            "COMBINED_INPUT_AUDIT_NOT_PASS",
            "Combined input audit status is not PASS",
            root_audit_path,
        )
        self.audit.expect(
            manifest.get("status") == "PASS",
            self.sample,
            stage,
            "INPUT_MANIFEST_NOT_PASS",
            "Sample input manifest status is not PASS",
            manifest_path,
        )
        combined_sample = combined.get("samples", {}).get(self.sample)
        self.audit.expect(
            isinstance(combined_sample, dict),
            self.sample,
            stage,
            "SAMPLE_ABSENT_FROM_COMBINED_AUDIT",
            "Sample is absent from combined input audit",
            root_audit_path,
        )
        if isinstance(combined_sample, dict):
            for key in ("records", "bp", "source_uncompressed_sha256", "short_fasta_sha256"):
                self.audit.expect(
                    combined_sample.get(key) == manifest.get(key),
                    self.sample,
                    stage,
                    "COMBINED_SAMPLE_MANIFEST_MISMATCH",
                    f"Combined audit and sample manifest disagree for {key}",
                    root_audit_path,
                    {"combined": combined_sample.get(key), "sample": manifest.get(key)},
                )

        expected_source = Path(self.config["source"])
        source = expected_source
        if not source.is_file() and Path(str(source) + ".gz").is_file():
            source = Path(str(source) + ".gz")
        if not self.require(source, stage, "locked source assembly"):
            return

        manifest_source = str(manifest.get("source_path", ""))
        allowed_source_strings = {str(expected_source), str(expected_source) + ".gz"}
        self.audit.expect(
            manifest_source in allowed_source_strings,
            self.sample,
            stage,
            "UNEXPECTED_SOURCE_PATH",
            "Input manifest source path is not the locked assembly path",
            manifest_path,
            {"observed": manifest_source, "allowed": sorted(allowed_source_strings)},
        )

        source_hash = sha256_uncompressed(source)
        self.audit.expect(
            source_hash == self.config["sha256"],
            self.sample,
            stage,
            "LOCKED_SOURCE_SHA256_MISMATCH",
            "Locked source assembly SHA-256 does not match",
            source,
            {"observed": source_hash, "expected": self.config["sha256"]},
        )
        self.audit.expect(
            manifest.get("source_uncompressed_sha256") == self.config["sha256"],
            self.sample,
            stage,
            "MANIFEST_SOURCE_SHA256_MISMATCH",
            "Input manifest does not contain the locked source SHA-256",
            manifest_path,
        )

        short_hash = sha256_file(short_path)
        self.audit.expect(
            short_hash == manifest.get("short_fasta_sha256"),
            self.sample,
            stage,
            "SHORT_FASTA_SHA256_MISMATCH",
            "Short-ID FASTA SHA-256 differs from its manifest",
            short_path,
            {"observed": short_hash, "manifest": manifest.get("short_fasta_sha256")},
        )

        with map_path.open("rt", encoding="utf-8", newline="") as handle:
            mapping_rows = list(csv.DictReader(handle, delimiter="\t"))
        required_columns = {"short_id", "original_id", "length", "sequence_sha256"}
        observed_columns = set(mapping_rows[0].keys()) if mapping_rows else set()
        self.audit.expect(
            required_columns <= observed_columns,
            self.sample,
            stage,
            "INVALID_CONTIG_MAP_HEADER",
            "Contig map lacks required columns",
            map_path,
            {"observed": sorted(observed_columns), "required": sorted(required_columns)},
        )

        count = 0
        total_bp = 0
        short_seen: Set[str] = set()
        original_seen: Set[str] = set()
        sequence_digest = hashlib.sha256()
        map_count = len(mapping_rows)
        for index, pair in enumerate(
            itertools.zip_longest(iter_fasta(source), iter_fasta(short_path)), 1
        ):
            original, short = pair
            if original is None or short is None:
                self.audit.error(
                    self.sample,
                    stage,
                    "SOURCE_SHORT_RECORD_COUNT_MISMATCH",
                    "Locked source and short-ID FASTA have different record counts",
                    short_path,
                    {"first_unpaired_record_index": index},
                )
                break
            count += 1
            total_bp += len(original.sequence)
            if original.identifier in original_seen:
                self.audit.error(
                    self.sample,
                    stage,
                    "DUPLICATE_ORIGINAL_ID",
                    f"Duplicate source ID {original.identifier}",
                    source,
                )
            original_seen.add(original.identifier)
            if short.identifier in short_seen:
                self.audit.error(
                    self.sample,
                    stage,
                    "DUPLICATE_SHORT_ID",
                    f"Duplicate short ID {short.identifier}",
                    short_path,
                )
            short_seen.add(short.identifier)
            expected_short_id = f"{self.config['short_prefix']}{index:06d}"
            if short.identifier != expected_short_id:
                self.audit.error(
                    self.sample,
                    stage,
                    "SHORT_ID_ORDER_MISMATCH",
                    f"Record {index} has short ID {short.identifier}; expected {expected_short_id}",
                    short_path,
                )
            if not original.sequence or not short.sequence:
                self.audit.error(
                    self.sample,
                    stage,
                    "EMPTY_ASSEMBLY_SEQUENCE",
                    f"Empty sequence at record {index}",
                    short_path,
                )
            invalid_original = sorted(set(original.sequence.upper()) - DNA_ALPHABET)
            invalid_short = sorted(set(short.sequence.upper()) - DNA_ALPHABET)
            if invalid_original or invalid_short:
                self.audit.error(
                    self.sample,
                    stage,
                    "INVALID_ASSEMBLY_SYMBOL",
                    f"Unsupported IUPAC symbol at record {index}",
                    short_path,
                    {"source": invalid_original, "short": invalid_short},
                )
            if original.sequence.upper() != short.sequence.upper():
                self.audit.error(
                    self.sample,
                    stage,
                    "SOURCE_SHORT_SEQUENCE_MISMATCH",
                    f"Source and short-ID sequences differ at record {index}",
                    short_path,
                    {"source_id": original.identifier, "short_id": short.identifier},
                )
            self.short_order.append(short.identifier)
            self.short_lengths[short.identifier] = len(short.sequence)
            sequence_digest.update(short.identifier.encode("ascii"))
            sequence_digest.update(b"\0")
            sequence_digest.update(short.sequence.upper().encode("ascii"))
            sequence_digest.update(b"\0")

            if index <= map_count and required_columns <= observed_columns:
                row = mapping_rows[index - 1]
                expected_sequence_hash = hashlib.sha256(
                    original.sequence.upper().encode("ascii")
                ).hexdigest()
                expected_row = {
                    "short_id": short.identifier,
                    "original_id": original.identifier,
                    "length": str(len(original.sequence)),
                    "sequence_sha256": expected_sequence_hash,
                }
                for key, expected_value in expected_row.items():
                    if row.get(key) != expected_value:
                        self.audit.error(
                            self.sample,
                            stage,
                            "CONTIG_MAP_ROW_MISMATCH",
                            f"Contig map row {index} disagrees for {key}",
                            map_path,
                            {"observed": row.get(key), "expected": expected_value},
                        )

        self.audit.expect(
            map_count == count,
            self.sample,
            stage,
            "CONTIG_MAP_COUNT_MISMATCH",
            "Contig map row count differs from FASTA record count",
            map_path,
            {"map_rows": map_count, "fasta_records": count},
        )
        self.audit.expect(
            count == self.config["records"],
            self.sample,
            stage,
            "LOCKED_RECORD_COUNT_MISMATCH",
            "Assembly record count differs from the locked value",
            source,
            {"observed": count, "expected": self.config["records"]},
        )
        self.audit.expect(
            total_bp == self.config["bp"],
            self.sample,
            stage,
            "LOCKED_TOTAL_BP_MISMATCH",
            "Assembly total length differs from the locked value",
            source,
            {"observed": total_bp, "expected": self.config["bp"]},
        )
        for key, expected in (("records", count), ("bp", total_bp)):
            self.audit.expect(
                manifest.get(key) == expected,
                self.sample,
                stage,
                "INPUT_MANIFEST_COUNT_MISMATCH",
                f"Input manifest {key} does not match re-audited FASTA",
                manifest_path,
                {"observed": manifest.get(key), "expected": expected},
            )
        self.audit.expect(
            manifest.get("ordered_shortid_sequence_sha256") == sequence_digest.hexdigest(),
            self.sample,
            stage,
            "ORDERED_SEQUENCE_DIGEST_MISMATCH",
            "Ordered short-ID sequence digest differs from manifest",
            manifest_path,
        )
        self.metrics.update(
            input_records=count,
            input_bp=total_bp,
            locked_source_sha256=source_hash,
            short_fasta_sha256=short_hash,
        )

    def audit_repeatmasker(self) -> None:
        stage = "repeatmasker"
        input_path = self.base / "00_input" / f"{self.sample}_final.shortids.fasta"
        out = self.base / "02_repeatmasker"
        masked = out / f"{self.sample}_final.shortids.fasta.masked"
        rm_out = out / f"{self.sample}_final.shortids.fasta.out"
        gff = out / f"{self.sample}_final.shortids.fasta.out.gff"
        tbl_path = out / f"{self.sample}_final.shortids.fasta.tbl"
        if not all(
            (
                self.require(input_path, stage, "RepeatMasker input FASTA"),
                self.require(masked, stage, "soft-masked FASTA"),
                self.require(rm_out, stage, "RepeatMasker .out table"),
                self.require(gff, stage, "RepeatMasker GFF"),
                self.require(tbl_path, stage, "RepeatMasker summary table"),
            )
        ):
            return

        count = 0
        total_bp = 0
        input_lowercase = 0
        masked_lowercase = 0
        newly_lowercased = 0
        unexpectedly_uppercased = 0
        masked_order: List[str] = []
        for index, pair in enumerate(
            itertools.zip_longest(iter_fasta(input_path), iter_fasta(masked)), 1
        ):
            original, soft = pair
            if original is None or soft is None:
                self.audit.error(
                    self.sample,
                    stage,
                    "MASKED_RECORD_COUNT_MISMATCH",
                    "Input and masked FASTAs have different record counts",
                    masked,
                    {"first_unpaired_record_index": index},
                )
                break
            count += 1
            total_bp += len(soft.sequence)
            masked_order.append(soft.identifier)
            self.masked_lengths[soft.identifier] = len(soft.sequence)
            if original.identifier != soft.identifier:
                self.audit.error(
                    self.sample,
                    stage,
                    "MASKED_ID_ORDER_MISMATCH",
                    f"Input/masked record {index} IDs differ",
                    masked,
                    {"input": original.identifier, "masked": soft.identifier},
                )
            if len(original.sequence) != len(soft.sequence):
                self.audit.error(
                    self.sample,
                    stage,
                    "MASKED_LENGTH_MISMATCH",
                    f"Input/masked lengths differ for {original.identifier}",
                    masked,
                    {"input": len(original.sequence), "masked": len(soft.sequence)},
                )
                continue
            if original.sequence.upper() != soft.sequence.upper():
                self.audit.error(
                    self.sample,
                    stage,
                    "MASKED_SEQUENCE_IDENTITY_FAILURE",
                    f"Masked sequence differs beyond letter case for {original.identifier}",
                    masked,
                )
            invalid = sorted(set(soft.sequence.upper()) - DNA_ALPHABET)
            if invalid:
                self.audit.error(
                    self.sample,
                    stage,
                    "INVALID_MASKED_SYMBOL",
                    f"Unsupported symbols in masked record {soft.identifier}: {''.join(invalid)}",
                    masked,
                )
            for before, after in zip(original.sequence, soft.sequence):
                input_lowercase += int(before.isalpha() and before.islower())
                masked_lowercase += int(after.isalpha() and after.islower())
                newly_lowercased += int(
                    before.isalpha()
                    and before.isupper()
                    and after.isalpha()
                    and after.islower()
                )
                unexpectedly_uppercased += int(
                    before.isalpha()
                    and before.islower()
                    and after.isalpha()
                    and after.isupper()
                )
        self.audit.expect(
            masked_order == self.short_order,
            self.sample,
            stage,
            "MASKED_FASTA_ORDER_MISMATCH",
            "Masked FASTA does not preserve the short-ID FASTA record order",
            masked,
            {
                "input_records": len(self.short_order),
                "masked_records": len(masked_order),
            },
        )
        self.audit.expect(
            count == self.config["records"] and total_bp == self.config["bp"],
            self.sample,
            stage,
            "MASKED_FASTA_SIZE_MISMATCH",
            "Masked FASTA count/total length differs from locked input",
            masked,
            {
                "records": count,
                "bp": total_bp,
                "expected_records": self.config["records"],
                "expected_bp": self.config["bp"],
            },
        )
        self.audit.expect(
            unexpectedly_uppercased == 0,
            self.sample,
            stage,
            "REPEATMASKER_UNMASKED_INPUT_BASES",
            "RepeatMasker converted pre-existing lowercase sequence to uppercase",
            masked,
            {"positions": unexpectedly_uppercased},
        )

        rm_out_rows = 0
        with rm_out.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not re.match(r"^\s*\d+\s+", raw):
                    continue
                fields = raw.split()
                if len(fields) < 8:
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_OUT_FIELD_COUNT",
                        f"RepeatMasker .out data line {line_number} is too short",
                        rm_out,
                    )
                    continue
                seqid = fields[4]
                try:
                    start, end = int(fields[5]), int(fields[6])
                except ValueError:
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_OUT_NONINTEGER_COORDINATE",
                        f"RepeatMasker .out line {line_number} has invalid coordinates",
                        rm_out,
                    )
                    continue
                if seqid not in self.masked_lengths:
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_OUT_UNKNOWN_SEQID",
                        f"RepeatMasker .out line {line_number} uses unknown seqid {seqid}",
                        rm_out,
                    )
                    continue
                if not (1 <= start <= end <= self.masked_lengths[seqid]):
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_OUT_OUT_OF_BOUNDS",
                        f"RepeatMasker .out interval on line {line_number} is out of bounds",
                        rm_out,
                        {
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                            "length": self.masked_lengths[seqid],
                        },
                    )
                    continue
                rm_out_rows += 1
        self.audit.expect(
            rm_out_rows > 0,
            self.sample,
            stage,
            "NO_REPEATMASKER_OUT_ROWS",
            "RepeatMasker .out contains no parseable repeat records",
            rm_out,
        )

        interval_by_contig: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        gff_features = 0
        with gff.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_GFF_FIELD_COUNT",
                        f"GFF line {line_number} has {len(fields)} fields, expected 9",
                        gff,
                    )
                    continue
                seqid = fields[0]
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError:
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_GFF_NONINTEGER_COORDINATE",
                        f"Non-integer coordinate on GFF line {line_number}",
                        gff,
                    )
                    continue
                if seqid not in self.masked_lengths:
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_GFF_UNKNOWN_SEQID",
                        f"Unknown masked FASTA seqid {seqid} on GFF line {line_number}",
                        gff,
                    )
                    continue
                if not (1 <= start <= end <= self.masked_lengths[seqid]):
                    self.audit.error(
                        self.sample,
                        stage,
                        "REPEATMASKER_GFF_OUT_OF_BOUNDS",
                        f"Out-of-bounds interval on GFF line {line_number}",
                        gff,
                        {
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                            "length": self.masked_lengths[seqid],
                        },
                    )
                    continue
                interval_by_contig[seqid].append((start, end))
                gff_features += 1
        self.audit.expect(
            gff_features > 0,
            self.sample,
            stage,
            "NO_REPEATMASKER_GFF_FEATURES",
            "RepeatMasker GFF contains no data features",
            gff,
        )

        union_bp = 0
        for intervals in interval_by_contig.values():
            current_start = current_end = None
            for start, end in sorted(intervals):
                if current_start is None:
                    current_start, current_end = start, end
                elif start <= current_end + 1:
                    current_end = max(current_end, end)
                else:
                    union_bp += current_end - current_start + 1
                    current_start, current_end = start, end
            if current_start is not None:
                union_bp += current_end - current_start + 1

        tbl_masked_bp: Optional[int] = None
        tbl_masked_percent: Optional[float] = None
        text = tbl_path.read_text(encoding="utf-8", errors="replace")
        match = re.search(
            r"bases masked:\s*([0-9,]+)\s*bp\s*\(\s*([0-9.]+)\s*%\)",
            text,
            re.IGNORECASE,
        )
        if match:
            tbl_masked_bp = int(match.group(1).replace(",", ""))
            tbl_masked_percent = float(match.group(2))
            if input_lowercase == 0:
                self.audit.expect(
                    tbl_masked_bp == masked_lowercase,
                    self.sample,
                    stage,
                    "REPEATMASKER_TBL_MASKED_BP_MISMATCH",
                    "RepeatMasker .tbl masked bp differs from lowercase bp in masked FASTA",
                    tbl_path,
                    {
                        "tbl_masked_bp": tbl_masked_bp,
                        "masked_fasta_lowercase_bp": masked_lowercase,
                    },
                )
        else:
            self.audit.error(
                self.sample,
                stage,
                "REPEATMASKER_TBL_UNPARSED",
                "RepeatMasker .tbl masked-bp line could not be parsed",
                tbl_path,
            )

        self.metrics.update(
            repeatmasker_input_records=len(self.short_order),
            repeatmasker_input_bp=sum(self.short_lengths.values()),
            repeatmasker_masked_records=count,
            repeatmasker_masked_bp=total_bp,
            repeatmasker_input_lowercase_bp=input_lowercase,
            repeatmasker_masked_lowercase_bp=masked_lowercase,
            repeatmasker_newly_lowercased_bp=newly_lowercased,
            repeatmasker_masked_percent=(100.0 * masked_lowercase / total_bp if total_bp else 0.0),
            repeatmasker_out_rows=rm_out_rows,
            repeatmasker_gff_features=gff_features,
            repeatmasker_gff_interval_union_bp=union_bp,
            repeatmasker_tbl_masked_bp=tbl_masked_bp,
            repeatmasker_tbl_masked_percent=tbl_masked_percent,
        )

    def audit_genome_busco(self) -> None:
        run_root = self.base / "01_genome_busco" / f"{self.sample}_final_genome"
        row = audit_busco_run(
            self.audit,
            self.sample,
            "genome_busco",
            run_root,
            expected_mode="genome",
            expected_input=self.base
            / "00_input"
            / f"{self.sample}_final.shortids.fasta",
            allowed_sequence_ids=set(self.short_lengths),
        )
        if row:
            self.audit.busco_rows.append(row)
            for key in (
                "complete", "single", "duplicated", "fragmented", "missing", "total",
                "complete_percent", "single_percent", "duplicated_percent",
                "fragmented_percent", "missing_percent",
            ):
                self.metrics[f"genome_busco_{key}"] = row[key]

    def audit_braker(self) -> None:
        stage = "braker"
        run = self.base / "03_braker4"
        samples_csv = run / "samples.csv"
        config_ini = run / "config.ini"
        dryrun_log = run / "dryrun.log"
        result = run / "output" / f"{self.sample}_final" / "results"
        gff = result / "braker.gff3.gz"
        proteins_path = result / "braker.aa.gz"
        cds_path = result / "braker.codingseq.gz"
        support_path = result / "gene_support.tsv"
        gene_stats_path = result / "quality_control" / "gene_set_statistics.txt"
        software_versions_path = result / "software_versions.tsv"
        citations_path = result / "braker_citations.bib"
        genome_path = result / "genome_path.txt"
        required = [
            self.require(samples_csv, stage, "BRAKER samples.csv"),
            self.require(config_ini, stage, "BRAKER config.ini"),
            self.require(dryrun_log, stage, "BRAKER Snakemake dry-run log"),
            self.require(gff, stage, "BRAKER GFF3"),
            self.require(proteins_path, stage, "BRAKER protein FASTA"),
            self.require(cds_path, stage, "BRAKER coding-sequence FASTA"),
            self.require(support_path, stage, "BRAKER gene-support table"),
            self.require(gene_stats_path, stage, "BRAKER gene-set statistics"),
            self.require(software_versions_path, stage, "BRAKER software versions"),
            self.require(citations_path, stage, "BRAKER citations"),
            self.require(genome_path, stage, "BRAKER genome provenance pointer"),
        ]
        if not all(required):
            return
        for path in (gff, proteins_path, cds_path):
            try:
                gzip_integrity(path)
            except Exception as exc:
                self.audit.error(
                    self.sample,
                    stage,
                    "BRAKER_GZIP_INTEGRITY_FAILURE",
                    f"gzip integrity check failed: {exc}",
                    path,
                )
                return

        with samples_csv.open("rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != 1:
            self.audit.error(
                self.sample,
                stage,
                "BRAKER_SAMPLE_ROW_COUNT",
                "BRAKER samples.csv must contain exactly one sample row",
                samples_csv,
                {"rows": len(rows)},
            )
        else:
            row = rows[0]
            expected_values = {
                "sample_name": f"{self.sample}_final",
                "genome": str(
                    self.base / "00_input" / f"{self.sample}_final.shortids.fasta"
                ),
                "genome_masked": str(
                    self.base
                    / "02_repeatmasker"
                    / f"{self.sample}_final.shortids.fasta.masked"
                ),
                "busco_lineage": EXPECTED_BUSCO_LINEAGE,
            }
            for key, expected in expected_values.items():
                self.audit.expect(
                    row.get(key) == expected,
                    self.sample,
                    stage,
                    "BRAKER_SAMPLE_CONFIG_MISMATCH",
                    f"BRAKER samples.csv has unexpected {key}",
                    samples_csv,
                    {"observed": row.get(key), "expected": expected},
                )
            rna_fields = (
                "bam_files",
                "fastq_r1",
                "fastq_r2",
                "sra_ids",
                "isoseq_bam",
                "isoseq_fastq",
            )
            populated_rna = {key: row.get(key) for key in rna_fields if row.get(key, "").strip()}
            self.audit.expect(
                not populated_rna,
                self.sample,
                stage,
                "UNEXPECTED_BRAKER_RNA_INPUT",
                "RNA evidence was not authorized but samples.csv contains RNA input",
                samples_csv,
                populated_rna,
            )
            self.metrics["braker_protein_evidence"] = row.get("protein_fasta", "")

        config_text = config_ini.read_text(encoding="utf-8")
        config_expectations = {
            "min_contig": str(self.config["braker_min_contig"]),
            "use_compleasm_hints": "1",
            "run_best_by_compleasm": "1",
            "skip_busco": "1",
            "translation_table": "1",
            "masking_tool": "repeatmasker",
        }
        for key, expected in config_expectations.items():
            match = re.search(
                rf"^{re.escape(key)}\s*=\s*(\S+)\s*$", config_text, re.MULTILINE
            )
            self.audit.expect(
                bool(match) and match.group(1) == expected,
                self.sample,
                stage,
                "BRAKER_CONFIG_MISMATCH",
                f"BRAKER config {key} is absent or differs from the submitted value",
                config_ini,
                {"observed": match.group(1) if match else None, "expected": expected},
            )

        genes: Dict[str, Tuple[str, int, int, str]] = {}
        transcripts: Dict[str, Tuple[str, str, int, int, str]] = {}
        cds_parents: Set[str] = set()
        gene_parent_ids: Set[str] = set()
        feature_counts: Counter = Counter()
        reversed_introns: List[dict] = []
        with gzip.open(gff, "rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_GFF_FIELD_COUNT",
                        f"GFF3 line {line_number} has {len(fields)} fields, expected 9",
                        gff,
                    )
                    continue
                seqid, _, feature, start_text, end_text, _, strand, phase, attr_text = fields
                feature_counts[feature] += 1
                try:
                    start, end = int(start_text), int(end_text)
                except ValueError:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_GFF_NONINTEGER_COORDINATE",
                        f"Non-integer GFF3 coordinate on line {line_number}",
                        gff,
                    )
                    continue
                if seqid not in self.masked_lengths:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_GFF_UNKNOWN_SEQID",
                        f"GFF3 line {line_number} uses seqid absent from masked FASTA: {seqid}",
                        gff,
                    )
                elif (
                    feature == "intron"
                    and start > end
                    and 1 <= end <= self.masked_lengths[seqid]
                    and 1 <= start <= self.masked_lengths[seqid]
                ):
                    # Some BRAKER releases emit a redundant zero/negative-span
                    # intron between abutting or overlapping exons.  It is not
                    # part of the submitted gene model (gene/mRNA/exon/CDS
                    # coordinates remain authoritative), so retain it as a
                    # visible warning rather than rejecting an otherwise valid
                    # result package.
                    reversed_introns.append(
                        {
                            "line": line_number,
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                        }
                    )
                elif not (1 <= start <= end <= self.masked_lengths[seqid]):
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_GFF_OUT_OF_BOUNDS",
                        f"GFF3 interval on line {line_number} is out of bounds",
                        gff,
                        {
                            "seqid": seqid,
                            "start": start,
                            "end": end,
                            "length": self.masked_lengths[seqid],
                        },
                    )
                if strand not in {"+", "-", ".", "?"}:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_GFF_INVALID_STRAND",
                        f"Invalid strand {strand!r} on line {line_number}",
                        gff,
                    )
                if feature == "CDS" and phase not in {"0", "1", "2"}:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_GFF_INVALID_CDS_PHASE",
                        f"Invalid CDS phase {phase!r} on line {line_number}",
                        gff,
                    )
                attrs = {}
                for part in attr_text.split(";"):
                    if "=" in part:
                        key, value = part.split("=", 1)
                        attrs[unquote(key)] = unquote(value)
                if feature == "gene":
                    gene_id = attrs.get("ID", "")
                    if not gene_id:
                        self.audit.error(
                            self.sample,
                            stage,
                            "BRAKER_GFF_GENE_WITHOUT_ID",
                            f"Gene on line {line_number} has no ID",
                            gff,
                        )
                    elif gene_id in genes:
                        self.audit.error(
                            self.sample,
                            stage,
                            "BRAKER_GFF_DUPLICATE_GENE_ID",
                            f"Duplicate gene ID {gene_id}",
                            gff,
                        )
                    else:
                        genes[gene_id] = (seqid, start, end, strand)
                elif feature in {"mRNA", "transcript"}:
                    transcript_id = attrs.get("ID", "")
                    parent_values = [x for x in attrs.get("Parent", "").split(",") if x]
                    if not transcript_id or len(parent_values) != 1:
                        self.audit.error(
                            self.sample,
                            stage,
                            "BRAKER_GFF_INVALID_TRANSCRIPT_RELATION",
                            f"Transcript on line {line_number} needs one ID and one gene Parent",
                            gff,
                            {"ID": transcript_id, "Parent": parent_values},
                        )
                    elif transcript_id in transcripts:
                        self.audit.error(
                            self.sample,
                            stage,
                            "BRAKER_GFF_DUPLICATE_TRANSCRIPT_ID",
                            f"Duplicate transcript ID {transcript_id}",
                            gff,
                        )
                    else:
                        transcripts[transcript_id] = (
                            parent_values[0],
                            seqid,
                            start,
                            end,
                            strand,
                        )
                        gene_parent_ids.add(parent_values[0])
                elif feature == "CDS":
                    parents = [x for x in attrs.get("Parent", "").split(",") if x]
                    if not parents:
                        self.audit.error(
                            self.sample,
                            stage,
                            "BRAKER_GFF_CDS_WITHOUT_PARENT",
                            f"CDS on line {line_number} has no Parent",
                            gff,
                        )
                    cds_parents.update(parents)

        if reversed_introns:
            self.audit.warning(
                self.sample,
                stage,
                "BRAKER_GFF_REVERSED_INTRON",
                "BRAKER GFF3 contains reversed redundant intron intervals; gene, transcript, exon and CDS coordinates are still audited strictly",
                gff,
                {
                    "count": len(reversed_introns),
                    "examples": reversed_introns[:20],
                },
            )

        self.audit.expect(
            bool(genes) and bool(transcripts) and feature_counts.get("CDS", 0) > 0,
            self.sample,
            stage,
            "EMPTY_BRAKER_MODEL_SET",
            "BRAKER GFF3 must contain genes, transcripts and CDS features",
            gff,
            dict(feature_counts),
        )
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "BRAKER_UNKNOWN_GENE_PARENT",
            "transcript_gene_parents",
            gene_parent_ids,
            "gene_ids",
            set(genes),
            gff,
        )
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "BRAKER_TRANSCRIPTS_WITHOUT_CDS",
            "CDS_parent_transcripts",
            cds_parents,
            "transcript_ids",
            set(transcripts),
            gff,
        )

        support_header: Optional[List[str]] = None
        support_rows: Dict[str, dict] = {}
        with support_path.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if support_header is None:
                    support_header = fields
                    continue
                if len(fields) != len(support_header):
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_SUPPORT_ROW_WIDTH_MISMATCH",
                        f"Gene-support line {line_number} has {len(fields)} fields; expected {len(support_header)}",
                        support_path,
                    )
                    continue
                row = dict(zip(support_header, fields))
                transcript_id = row.get("transcript_id", "")
                if not transcript_id:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_SUPPORT_MISSING_TRANSCRIPT_ID",
                        f"Gene-support line {line_number} has no transcript_id",
                        support_path,
                    )
                elif transcript_id in support_rows:
                    self.audit.error(
                        self.sample,
                        stage,
                        "BRAKER_SUPPORT_DUPLICATE_TRANSCRIPT_ID",
                        f"Duplicate gene-support transcript {transcript_id}",
                        support_path,
                    )
                else:
                    support_rows[transcript_id] = row
        required_support_columns = {"gene_id", "transcript_id", "chrom", "strand"}
        self.audit.expect(
            support_header is not None
            and required_support_columns <= set(support_header),
            self.sample,
            stage,
            "BRAKER_SUPPORT_HEADER_INVALID",
            "Gene-support table lacks required structural columns",
            support_path,
            {
                "observed": support_header,
                "required": sorted(required_support_columns),
            },
        )
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "BRAKER_SUPPORT_TRANSCRIPT_SET_MISMATCH",
            "support_transcript_ids",
            set(support_rows),
            "GFF_transcript_ids",
            set(transcripts),
            support_path,
        )
        support_structure_mismatch: List[str] = []
        for transcript_id in set(support_rows) & set(transcripts):
            gene_id, seqid, _, _, strand = transcripts[transcript_id]
            row = support_rows[transcript_id]
            if (
                row.get("gene_id") != gene_id
                or row.get("chrom") != seqid
                or row.get("strand") != strand
            ):
                support_structure_mismatch.append(transcript_id)
        if support_structure_mismatch:
            self.audit.error(
                self.sample,
                stage,
                "BRAKER_SUPPORT_STRUCTURE_MISMATCH",
                "Gene-support gene/contig/strand fields disagree with GFF3",
                support_path,
                {
                    "count": len(support_structure_mismatch),
                    "examples": sorted(support_structure_mismatch)[:20],
                },
            )

        proteins, protein_order = read_fasta_map(
            proteins_path,
            self.audit,
            self.sample,
            stage,
            PROTEIN_ALPHABET,
            "protein",
        )
        cds, _ = read_fasta_map(
            cds_path,
            self.audit,
            self.sample,
            stage,
            DNA_ALPHABET,
            "coding DNA",
        )
        tx_ids = set(transcripts)
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "BRAKER_GFF_PROTEIN_ID_MISMATCH",
            "GFF_transcripts",
            tx_ids,
            "protein_ids",
            set(proteins),
            proteins_path,
        )
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "BRAKER_GFF_CDS_ID_MISMATCH",
            "GFF_transcripts",
            tx_ids,
            "codingseq_ids",
            set(cds),
            cds_path,
        )
        internal_stop_ids = [
            identifier
            for identifier, sequence in proteins.items()
            if "*" in sequence.rstrip("*")
        ]
        if internal_stop_ids:
            self.audit.error(
                self.sample,
                stage,
                "BRAKER_PROTEIN_INTERNAL_STOPS",
                "Predicted proteins contain internal stop codons",
                proteins_path,
                {"count": len(internal_stop_ids), "examples": internal_stop_ids[:20]},
            )

        translation_mismatches: List[str] = []
        cds_modulo_errors: List[str] = []
        for identifier in sorted(set(proteins) & set(cds)):
            cds_sequence = cds[identifier]
            if len(cds_sequence) % 3:
                cds_modulo_errors.append(identifier)
            translated = translate_dna(cds_sequence).rstrip("*")
            if translated != proteins[identifier].rstrip("*"):
                translation_mismatches.append(identifier)
        if cds_modulo_errors:
            self.audit.error(
                self.sample,
                stage,
                "BRAKER_CDS_NOT_MODULO_THREE",
                "BRAKER coding sequences are not all divisible by three",
                cds_path,
                {"count": len(cds_modulo_errors), "examples": cds_modulo_errors[:20]},
            )
        if translation_mismatches:
            self.audit.error(
                self.sample,
                stage,
                "BRAKER_CDS_TRANSLATION_MISMATCH",
                "Coding-sequence translations do not match BRAKER proteins",
                cds_path,
                {
                    "count": len(translation_mismatches),
                    "examples": translation_mismatches[:20],
                },
            )

        self.gff_gene_ids = set(genes)
        self.gff_transcript_ids = set(transcripts)
        self.tx_to_gene = {key: value[0] for key, value in transcripts.items()}
        self.braker_proteins = proteins
        self.braker_protein_order = protein_order
        self.braker_cds = cds
        self.metrics.update(
            braker_genes=len(genes),
            braker_transcripts=len(transcripts),
            braker_cds_features=feature_counts.get("CDS", 0),
            braker_protein_isoforms=len(proteins),
            braker_coding_sequences=len(cds),
            braker_internal_stop_proteins=len(internal_stop_ids),
            braker_cds_translation_mismatches=len(translation_mismatches),
            braker_gene_support_rows=len(support_rows),
            braker_gene_support_structure_mismatches=len(
                support_structure_mismatch
            ),
            braker_reversed_intron_features=len(reversed_introns),
            braker_gff_feature_counts=dict(feature_counts),
            braker_uses_compleasm_hints=True,
            braker_has_declared_rna_evidence=False,
        )
        for optional_name in ("training_summary.txt", "compleasm_summary.txt"):
            optional_path = result / optional_name
            if optional_path.is_file() and optional_path.stat().st_size > 0:
                self.audit.add_file(
                    optional_path, self.sample, stage, f"optional BRAKER {optional_name}"
                )

    def audit_protein_and_busco(self) -> None:
        stage = "protein_busco"
        out = self.base / "04_protein_busco"
        all_path = out / f"{self.sample}_final.all_isoforms.aa"
        representative_path = out / f"{self.sample}_final.longest_isoform.aa"
        stats_path = out / "isoform_stats.tsv"
        complete_marker = out / ".complete"
        braker_gz = (
            self.base
            / "03_braker4"
            / "output"
            / f"{self.sample}_final"
            / "results"
            / "braker.aa.gz"
        )
        if not all(
            (
                self.require(all_path, stage, "all-isoform protein FASTA"),
                self.require(representative_path, stage, "representative protein FASTA"),
                self.require(stats_path, stage, "isoform statistics"),
                self.require(complete_marker, stage, "protein BUSCO completion marker", nonempty=False),
                self.require(braker_gz, stage, "BRAKER protein FASTA"),
            )
        ):
            return

        decompressed_hash = sha256_uncompressed(braker_gz)
        all_hash = sha256_file(all_path)
        self.audit.expect(
            decompressed_hash == all_hash,
            self.sample,
            stage,
            "ALL_ISOFORMS_NOT_EXACT_BRAKER_COPY",
            "All-isoform FASTA is not the exact decompressed BRAKER protein file",
            all_path,
            {"braker_uncompressed_sha256": decompressed_hash, "all_sha256": all_hash},
        )
        all_proteins, all_order = read_fasta_map(
            all_path,
            self.audit,
            self.sample,
            stage,
            PROTEIN_ALPHABET,
            "all-isoform protein",
        )
        representatives, representative_order = read_fasta_map(
            representative_path,
            self.audit,
            self.sample,
            stage,
            PROTEIN_ALPHABET,
            "representative protein",
        )
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "ALL_ISOFORMS_BRAKER_ID_MISMATCH",
            "all_isoform_ids",
            set(all_proteins),
            "BRAKER_protein_ids",
            set(self.braker_proteins),
            all_path,
        )
        sequence_mismatch = [
            identifier
            for identifier in set(all_proteins) & set(self.braker_proteins)
            if all_proteins[identifier] != self.braker_proteins[identifier]
        ]
        if sequence_mismatch:
            self.audit.error(
                self.sample,
                stage,
                "ALL_ISOFORMS_BRAKER_SEQUENCE_MISMATCH",
                "All-isoform sequences differ from BRAKER proteins",
                all_path,
                {"count": len(sequence_mismatch), "examples": sorted(sequence_mismatch)[:20]},
            )
        self.audit.expect(
            0 < len(representatives) <= len(all_proteins),
            self.sample,
            stage,
            "INVALID_REPRESENTATIVE_COUNT",
            "Representative protein count must be positive and no larger than all isoforms",
            representative_path,
            {"representatives": len(representatives), "all_isoforms": len(all_proteins)},
        )
        unknown_transcripts = [identifier for identifier in all_order if identifier not in self.tx_to_gene]
        if unknown_transcripts:
            self.audit.error(
                self.sample,
                stage,
                "ISOFORMS_ABSENT_FROM_GFF",
                "All-isoform FASTA contains IDs absent from BRAKER GFF transcripts",
                all_path,
                {"count": len(unknown_transcripts), "examples": unknown_transcripts[:20]},
            )

        expected_ids: List[str] = []
        first_gene_order: List[str] = []
        best_by_gene: Dict[str, Tuple[str, int]] = {}
        for identifier in all_order:
            gene = self.tx_to_gene.get(identifier)
            if gene is None:
                continue
            length = len(all_proteins[identifier])
            if gene not in best_by_gene:
                first_gene_order.append(gene)
                best_by_gene[gene] = (identifier, length)
            elif length > best_by_gene[gene][1]:
                best_by_gene[gene] = (identifier, length)
        expected_ids = [best_by_gene[gene][0] for gene in first_gene_order]
        self.audit.expect(
            representative_order == expected_ids,
            self.sample,
            stage,
            "LONGEST_ISOFORM_SELECTION_MISMATCH",
            "Representative FASTA is not exactly one first-occurring longest isoform per GFF gene",
            representative_path,
            {
                "observed_count": len(representative_order),
                "expected_count": len(expected_ids),
                "first_observed": representative_order[:20],
                "first_expected": expected_ids[:20],
            },
        )
        representative_sequence_mismatch = [
            identifier
            for identifier in set(representatives) & set(all_proteins)
            if representatives[identifier] != all_proteins[identifier]
        ]
        if representative_sequence_mismatch:
            self.audit.error(
                self.sample,
                stage,
                "REPRESENTATIVE_SEQUENCE_MISMATCH",
                "Representative sequences differ from their all-isoform records",
                representative_path,
                {
                    "count": len(representative_sequence_mismatch),
                    "examples": sorted(representative_sequence_mismatch)[:20],
                },
            )
        self.audit.expect(
            len(representatives) == len(best_by_gene),
            self.sample,
            stage,
            "REPRESENTATIVE_GENE_COUNT_MISMATCH",
            "Representative count does not equal the number of genes represented in all isoforms",
            representative_path,
            {"representatives": len(representatives), "genes": len(best_by_gene)},
        )
        if self.gff_gene_ids:
            compare_sets(
                self.audit,
                self.sample,
                stage,
                "REPRESENTATIVE_GFF_GENE_MISMATCH",
                "represented_genes",
                set(best_by_gene),
                "GFF_gene_ids",
                self.gff_gene_ids,
                representative_path,
            )

        isoform_stats: Dict[str, int] = {}
        with stats_path.open("rt", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 2 or not fields[1].isdigit():
                    self.audit.error(
                        self.sample,
                        stage,
                        "INVALID_ISOFORM_STATS_ROW",
                        f"Invalid isoform_stats.tsv row {line_number}",
                        stats_path,
                    )
                    continue
                isoform_stats[fields[0]] = int(fields[1])
        self.audit.expect(
            isoform_stats.get("all_isoforms") == len(all_proteins),
            self.sample,
            stage,
            "ALL_ISOFORM_STATS_MISMATCH",
            "isoform_stats.tsv all_isoforms does not match FASTA",
            stats_path,
        )
        self.audit.expect(
            isoform_stats.get("representative_longest_isoforms") == len(representatives),
            self.sample,
            stage,
            "REPRESENTATIVE_STATS_MISMATCH",
            "isoform_stats.tsv representative count does not match FASTA",
            stats_path,
        )
        self.representatives = representatives
        self.representative_order = representative_order
        self.metrics.update(
            all_protein_isoforms=len(all_proteins),
            representative_longest_proteins=len(representatives),
            represented_genes=len(best_by_gene),
        )

        run_root = out / f"{self.sample}_final_longest_metazoa"
        row = audit_busco_run(
            self.audit,
            self.sample,
            stage,
            run_root,
            expected_mode="proteins",
            expected_input=representative_path,
            allowed_sequence_ids=set(representatives),
        )
        if row:
            self.audit.busco_rows.append(row)
            for key in (
                "complete", "single", "duplicated", "fragmented", "missing", "total",
                "complete_percent", "single_percent", "duplicated_percent",
                "fragmented_percent", "missing_percent",
            ):
                self.metrics[f"protein_busco_{key}"] = row[key]

    def audit_eggnog(self) -> None:
        stage = "eggnog"
        path = (
            self.base
            / "05_eggnog"
            / f"{self.sample}_final.longest.emapper.annotations"
        )
        if not self.require(path, stage, "eggNOG-mapper annotations"):
            return
        if not self.representatives:
            self.audit.error(
                self.sample,
                stage,
                "REPRESENTATIVE_INPUT_UNAVAILABLE",
                "Representative proteins were unavailable for eggNOG query-ID validation",
                path,
            )
            return
        header: Optional[List[str]] = None
        comments: List[str] = []
        query_ids: List[str] = []
        scanned_footer_values: List[int] = []
        seed_ortholog_count = 0
        row_widths: Counter = Counter()
        with path.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#query\t") or line == "#query":
                    candidate = line[1:].split("\t")
                    if header is not None and header != candidate:
                        self.audit.error(
                            self.sample,
                            stage,
                            "EGGNOG_INCONSISTENT_HEADER",
                            "eggNOG file contains multiple different #query headers",
                            path,
                        )
                    header = candidate
                    continue
                if line.startswith("#"):
                    comments.append(line)
                    footer_match = re.fullmatch(
                        r"##\s+(\d+)\s+queries scanned\s*", line
                    )
                    if footer_match:
                        scanned_footer_values.append(int(footer_match.group(1)))
                    continue
                fields = line.split("\t")
                row_widths[len(fields)] += 1
                if header is None:
                    self.audit.error(
                        self.sample,
                        stage,
                        "EGGNOG_DATA_BEFORE_HEADER",
                        f"eggNOG data row precedes #query header at line {line_number}",
                        path,
                    )
                    continue
                if len(fields) != len(header):
                    self.audit.error(
                        self.sample,
                        stage,
                        "EGGNOG_ROW_WIDTH_MISMATCH",
                        f"eggNOG line {line_number} has {len(fields)} fields; expected {len(header)}",
                        path,
                    )
                    continue
                row = dict(zip(header, fields))
                query = row.get("query", fields[0])
                query_ids.append(query)
                if row.get("seed_ortholog", "") not in MISSING_VALUES:
                    seed_ortholog_count += 1
        self.audit.expect(
            header is not None,
            self.sample,
            stage,
            "EGGNOG_HEADER_MISSING",
            "eggNOG annotations lack a #query header",
            path,
        )
        if header is not None:
            required = {"query", "seed_ortholog", "COG_category", "GOs", "KEGG_ko", "KEGG_Pathway"}
            self.audit.expect(
                required <= set(header),
                self.sample,
                stage,
                "EGGNOG_REQUIRED_COLUMNS_MISSING",
                "eggNOG header lacks required annotation columns",
                path,
                {"missing": sorted(required - set(header))},
            )
            self.audit.expect(
                len(header) == len(set(header)),
                self.sample,
                stage,
                "EGGNOG_DUPLICATE_COLUMNS",
                "eggNOG header contains duplicate column names",
                path,
            )
        self.audit.expect(
            bool(query_ids),
            self.sample,
            stage,
            "EGGNOG_NO_DATA_ROWS",
            "eggNOG annotations contain no data rows",
            path,
        )
        duplicate_queries = [key for key, value in Counter(query_ids).items() if value > 1]
        if duplicate_queries:
            self.audit.error(
                self.sample,
                stage,
                "EGGNOG_DUPLICATE_QUERY",
                "eggNOG annotations contain duplicate query rows",
                path,
                {"count": len(duplicate_queries), "examples": sorted(duplicate_queries)[:20]},
            )
        unknown = set(query_ids) - set(self.representatives)
        if unknown:
            self.audit.error(
                self.sample,
                stage,
                "EGGNOG_UNKNOWN_QUERY",
                "eggNOG query IDs are absent from representative proteins",
                path,
                {"count": len(unknown), "examples": sorted(unknown)[:20]},
            )
        self.audit.expect(
            len(scanned_footer_values) == 1,
            self.sample,
            stage,
            "EGGNOG_QUERY_SCAN_FOOTER_MISSING_OR_DUPLICATE",
            "eggNOG annotations must contain exactly one successful-run 'queries scanned' footer",
            path,
            {"observed_values": scanned_footer_values},
        )
        scanned_queries = (
            scanned_footer_values[0] if len(scanned_footer_values) == 1 else None
        )
        if scanned_queries is not None:
            self.audit.expect(
                len(set(query_ids)) <= scanned_queries <= len(self.representatives),
                self.sample,
                stage,
                "EGGNOG_QUERY_SCAN_COUNT_INCONSISTENT",
                "eggNOG scanned-query count must cover annotation rows and cannot exceed representative proteins",
                path,
                {
                    "annotation_queries": len(set(query_ids)),
                    "scanned_queries": scanned_queries,
                    "representative_proteins": len(self.representatives),
                },
            )
        comment_text = "\n".join(comments)
        version_match = re.search(r"emapper[-_ ](?:version[:= ]*)?v?(\d+\.\d+\.\d+)", comment_text, re.I)
        eggnog_version = version_match.group(1) if version_match else ""
        if eggnog_version:
            self.audit.expect(
                eggnog_version == "2.1.13",
                self.sample,
                stage,
                "EGGNOG_VERSION_MISMATCH",
                "eggNOG-mapper version differs from the pinned workflow",
                path,
                {"observed": eggnog_version, "expected": "2.1.13"},
            )
        else:
            self.audit.warning(
                self.sample,
                stage,
                "EGGNOG_VERSION_NOT_RECOVERED",
                "eggNOG-mapper version was not parsed from output comments",
                path,
            )
        if "--tax_scope 33208" not in comment_text and "--tax_scope=33208" not in comment_text:
            self.audit.warning(
                self.sample,
                stage,
                "EGGNOG_TAX_SCOPE_NOT_RECOVERED",
                "Metazoa tax_scope 33208 was not found in output comments",
                path,
            )
        self.metrics.update(
            eggnog_annotation_rows=len(query_ids),
            eggnog_unique_queries=len(set(query_ids)),
            eggnog_queries_scanned=scanned_queries,
            eggnog_seed_ortholog_queries=seed_ortholog_count,
            eggnog_seed_ortholog_percent=(
                100.0 * seed_ortholog_count / len(self.representatives)
                if self.representatives
                else 0.0
            ),
            eggnog_mapper_version=eggnog_version,
        )
        self.audit.functional_rows.append(
            {
                "sample": self.sample,
                "source": "eggNOG",
                "data_rows": len(query_ids),
                "unique_queries": len(set(query_ids)),
                "representative_proteins": len(self.representatives),
                "queries_scanned": scanned_queries,
                "annotated_queries": seed_ortholog_count,
                "annotated_percent": (
                    100.0 * seed_ortholog_count / len(self.representatives)
                    if self.representatives
                    else 0.0
                ),
                "details": {"row_widths": dict(row_widths)},
            }
        )

    def audit_swissprot(self) -> None:
        stage = "swissprot"
        out = self.base / "06_swissprot"
        all_path = out / f"{self.sample}_final_vs_swissprot.tsv.gz"
        best_path = out / f"{self.sample}_final_vs_swissprot.best.tsv"
        if not all(
            (
                self.require(all_path, stage, "all Swiss-Prot DIAMOND hits"),
                self.require(best_path, stage, "Swiss-Prot best hits"),
            )
        ):
            return
        try:
            gzip_integrity(all_path)
        except Exception as exc:
            self.audit.error(
                self.sample,
                stage,
                "SWISSPROT_GZIP_INTEGRITY_FAILURE",
                f"Swiss-Prot all-hit gzip integrity failed: {exc}",
                all_path,
            )
            return
        if not self.representatives:
            self.audit.error(
                self.sample,
                stage,
                "REPRESENTATIVE_INPUT_UNAVAILABLE",
                "Representative proteins were unavailable for Swiss-Prot validation",
                all_path,
            )
            return
        all_rows: List[List[str]] = []
        by_query: Dict[str, List[List[str]]] = defaultdict(list)
        with gzip.open(all_path, "rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 11:
                    self.audit.error(
                        self.sample,
                        stage,
                        "SWISSPROT_FIELD_COUNT",
                        f"All-hit line {line_number} has {len(fields)} fields, expected 11",
                        all_path,
                    )
                    continue
                validate_swissprot_fields(
                    self.audit,
                    self.sample,
                    stage,
                    all_path,
                    line_number,
                    fields,
                    self.representatives,
                )
                all_rows.append(fields)
                by_query[fields[0]].append(fields)
        best_rows: List[List[str]] = []
        with best_path.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 11:
                    self.audit.error(
                        self.sample,
                        stage,
                        "SWISSPROT_BEST_FIELD_COUNT",
                        f"Best-hit line {line_number} has {len(fields)} fields, expected 11",
                        best_path,
                    )
                    continue
                validate_swissprot_fields(
                    self.audit,
                    self.sample,
                    stage,
                    best_path,
                    line_number,
                    fields,
                    self.representatives,
                )
                best_rows.append(fields)
        self.audit.expect(
            bool(all_rows) and bool(best_rows),
            self.sample,
            stage,
            "SWISSPROT_NO_DATA_ROWS",
            "Swiss-Prot all-hit and best-hit tables must contain data rows",
            out,
        )
        too_many = {
            query: len({row[1] for row in rows})
            for query, rows in by_query.items()
            if len({row[1] for row in rows}) > 5
        }
        if too_many:
            self.audit.error(
                self.sample,
                stage,
                "SWISSPROT_MAX_TARGETS_EXCEEDED",
                "Swiss-Prot output has more than five targets for a query",
                all_path,
                dict(list(sorted(too_many.items()))[:20]),
            )
        best_by_query: Dict[str, List[str]] = {}
        for fields in best_rows:
            if fields[0] in best_by_query:
                self.audit.error(
                    self.sample,
                    stage,
                    "SWISSPROT_DUPLICATE_BEST_QUERY",
                    f"Best-hit table contains multiple rows for {fields[0]}",
                    best_path,
                )
            else:
                best_by_query[fields[0]] = fields
        compare_sets(
            self.audit,
            self.sample,
            stage,
            "SWISSPROT_BEST_QUERY_SET_MISMATCH",
            "best_query_ids",
            set(best_by_query),
            "all_hit_query_ids",
            set(by_query),
            best_path,
        )
        wrong_first: List[str] = []
        first_not_max_reported_bitscore = 0
        for query, rows in by_query.items():
            if query not in best_by_query:
                continue
            if best_by_query[query] != rows[0]:
                wrong_first.append(query)
            try:
                maximum_reported_bitscore = max(float(row[9]) for row in rows)
                if not math.isclose(
                    float(rows[0][9]), maximum_reported_bitscore
                ):
                    first_not_max_reported_bitscore += 1
            except ValueError:
                pass
        if wrong_first:
            self.audit.error(
                self.sample,
                stage,
                "SWISSPROT_BEST_NOT_FIRST",
                "Best-hit rows are not exact first rows from the all-hit table",
                best_path,
                {"count": len(wrong_first), "examples": wrong_first[:20]},
            )
        self.metrics.update(
            swissprot_all_hit_rows=len(all_rows),
            swissprot_best_hit_queries=len(best_by_query),
            swissprot_first_not_max_reported_bitscore_queries=(
                first_not_max_reported_bitscore
            ),
            swissprot_best_hit_percent=(
                100.0 * len(best_by_query) / len(self.representatives)
                if self.representatives
                else 0.0
            ),
        )
        self.audit.functional_rows.append(
            {
                "sample": self.sample,
                "source": "Swiss-Prot",
                "data_rows": len(all_rows),
                "unique_queries": len(by_query),
                "representative_proteins": len(self.representatives),
                "annotated_queries": len(best_by_query),
                "annotated_percent": (
                    100.0 * len(best_by_query) / len(self.representatives)
                    if self.representatives
                    else 0.0
                ),
                "details": {"maximum_hits_per_query": max(map(len, by_query.values()), default=0)},
            }
        )

    def audit_interpro(self) -> None:
        stage = "interpro"
        out = self.base / "07_interpro"
        plain = out / f"{self.sample}_final.longest.interproscan.tsv"
        compressed = Path(str(plain) + ".gz")
        if not compressed.is_file():
            self.audit.error(
                self.sample,
                stage,
                "MISSING_FILE",
                "Missing required compressed InterProScan TSV.gz final output",
                compressed,
                {"uncompressed_output_exists": plain.is_file()},
            )
            return
        if plain.is_file():
            self.audit.error(
                self.sample,
                stage,
                "AMBIGUOUS_INTERPRO_FINAL_OUTPUT",
                "Both uncompressed and compressed InterProScan final outputs exist",
                out,
                [str(plain), str(compressed)],
            )
        path = compressed
        if not self.require(path, stage, "InterProScan TSV"):
            return
        try:
            gzip_integrity(path)
        except Exception as exc:
            self.audit.error(
                self.sample,
                stage,
                "INTERPRO_GZIP_INTEGRITY_FAILURE",
                f"InterProScan gzip integrity failed: {exc}",
                path,
            )
            return
        if not self.representatives:
            self.audit.error(
                self.sample,
                stage,
                "REPRESENTATIVE_INPUT_UNAVAILABLE",
                "Representative proteins were unavailable for InterProScan validation",
                path,
            )
            return
        query_ids: Set[str] = set()
        ipr_query_ids: Set[str] = set()
        analyses: Set[str] = set()
        row_count = 0
        row_widths: Counter = Counter()
        with open_text(path) as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                row_count += 1
                row_widths[len(fields)] += 1
                if not 11 <= len(fields) <= 15:
                    self.audit.error(
                        self.sample,
                        stage,
                        "INTERPRO_FIELD_COUNT",
                        f"InterProScan line {line_number} has {len(fields)} fields; expected 11-15",
                        path,
                    )
                    continue
                query = fields[0]
                query_ids.add(query)
                if query not in self.representatives:
                    self.audit.error(
                        self.sample,
                        stage,
                        "INTERPRO_UNKNOWN_QUERY",
                        f"InterProScan query {query} is absent from representative proteins",
                        path,
                    )
                    continue
                protein = self.representatives[query].rstrip("*")
                observed_md5 = fields[1].lower()
                expected_md5 = hashlib.md5(protein.encode("ascii")).hexdigest()
                if observed_md5 != expected_md5:
                    self.audit.error(
                        self.sample,
                        stage,
                        "INTERPRO_SEQUENCE_MD5_MISMATCH",
                        f"InterProScan sequence MD5 differs for {query}",
                        path,
                        {"observed": observed_md5, "expected": expected_md5},
                    )
                try:
                    reported_length = int(fields[2])
                    start, end = int(fields[6]), int(fields[7])
                except ValueError:
                    self.audit.error(
                        self.sample,
                        stage,
                        "INTERPRO_NONINTEGER_LENGTH_OR_COORDINATE",
                        f"Invalid length/domain coordinate on line {line_number}",
                        path,
                    )
                    continue
                if reported_length != len(protein):
                    self.audit.error(
                        self.sample,
                        stage,
                        "INTERPRO_SEQUENCE_LENGTH_MISMATCH",
                        f"InterProScan sequence length differs for {query}",
                        path,
                        {"observed": reported_length, "expected": len(protein)},
                    )
                if not (1 <= start <= end <= reported_length):
                    self.audit.error(
                        self.sample,
                        stage,
                        "INTERPRO_DOMAIN_OUT_OF_BOUNDS",
                        f"InterProScan domain coordinates are out of bounds on line {line_number}",
                        path,
                        {"start": start, "end": end, "length": reported_length},
                    )
                analyses.add(fields[3])
                if len(fields) >= 12 and fields[11] not in MISSING_VALUES:
                    if not IPR_RE.match(fields[11]):
                        self.audit.error(
                            self.sample,
                            stage,
                            "INTERPRO_INVALID_IPR_ACCESSION",
                            f"Invalid InterPro accession {fields[11]!r} on line {line_number}",
                            path,
                        )
                    else:
                        ipr_query_ids.add(query)
        self.audit.expect(
            row_count > 0,
            self.sample,
            stage,
            "INTERPRO_NO_DATA_ROWS",
            "InterProScan output contains no data rows",
            path,
        )
        self.metrics.update(
            interproscan_rows=row_count,
            interproscan_signature_queries=len(query_ids),
            interproscan_signature_percent=(
                100.0 * len(query_ids) / len(self.representatives)
                if self.representatives
                else 0.0
            ),
            interpro_annotated_queries=len(query_ids),
            interpro_annotated_percent=(
                100.0 * len(query_ids) / len(self.representatives)
                if self.representatives
                else 0.0
            ),
            interpro_integrated_ipr_queries=len(ipr_query_ids),
            interpro_integrated_ipr_percent=(
                100.0 * len(ipr_query_ids) / len(self.representatives)
                if self.representatives
                else 0.0
            ),
            interproscan_analyses=sorted(analyses),
        )
        self.audit.functional_rows.append(
            {
                "sample": self.sample,
                "source": "InterProScan",
                "data_rows": row_count,
                "unique_queries": len(query_ids),
                "representative_proteins": len(self.representatives),
                "annotated_queries": len(query_ids),
                "annotated_percent": (
                    100.0 * len(query_ids) / len(self.representatives)
                    if self.representatives
                    else 0.0
                ),
                "details": {
                    "integrated_interpro_queries": len(ipr_query_ids),
                    "integrated_interpro_percent": (
                        100.0 * len(ipr_query_ids) / len(self.representatives)
                        if self.representatives
                        else 0.0
                    ),
                    "row_widths": dict(row_widths),
                    "analyses": sorted(analyses),
                },
            }
        )


def parse_busco_summary(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="strict")
    patterns = {
        "complete": r"^\s*(\d+)\s+Complete BUSCOs",
        "single": r"^\s*(\d+)\s+Complete and single-copy BUSCOs",
        "duplicated": r"^\s*(\d+)\s+Complete and duplicated BUSCOs",
        "fragmented": r"^\s*(\d+)\s+Fragmented BUSCOs",
        "missing": r"^\s*(\d+)\s+Missing BUSCOs",
        "total": r"^\s*(\d+)\s+Total BUSCO groups searched",
    }
    result = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.MULTILINE)
        if not match:
            raise ValueError(f"Could not parse BUSCO {key} count from {path}")
        result[key] = int(match.group(1))
    version = re.search(r"BUSCO version is:\s*([^\s]+)", text)
    lineage = re.search(r"lineage dataset is:\s*([^\s(]+)", text, re.IGNORECASE)
    mode = re.search(r"BUSCO was run in mode:\s*([^\s]+)", text, re.IGNORECASE)
    input_match = re.search(
        r"Summarized benchmarking in BUSCO notation for file\s+(.+?)\s*$",
        text,
        re.MULTILINE,
    )
    result.update(
        version=version.group(1) if version else "",
        lineage=lineage.group(1) if lineage else "",
        mode=mode.group(1) if mode else "",
        input_path=input_match.group(1).strip() if input_match else "",
        path=str(path),
    )
    return result


def audit_busco_run(
    audit: Audit,
    sample: str,
    stage: str,
    run_root: Path,
    expected_mode: str,
    expected_input: Path,
    allowed_sequence_ids: Set[str],
) -> Optional[dict]:
    if not run_root.is_dir():
        audit.error(sample, stage, "MISSING_BUSCO_RUN", "Missing BUSCO run directory", run_root)
        return None
    summaries = sorted(
        set(run_root.glob("short_summary*.txt"))
        | set(
            (run_root / f"run_{EXPECTED_BUSCO_LINEAGE}").glob(
                "short_summary*.txt"
            )
        )
    )
    summaries = [path for path in summaries if path.is_file() and path.stat().st_size > 0]
    if not summaries:
        audit.error(sample, stage, "MISSING_BUSCO_SUMMARY", "No BUSCO short summary found", run_root)
        return None
    for path in summaries:
        audit.add_file(path, sample, stage, "BUSCO short summary")
    parsed = [parse_busco_summary(path) for path in summaries]
    canonical = parsed[0]
    count_keys = ("complete", "single", "duplicated", "fragmented", "missing", "total")
    for other in parsed[1:]:
        if any(other[key] != canonical[key] for key in count_keys):
            audit.error(
                sample,
                stage,
                "BUSCO_SUMMARIES_DISAGREE",
                "BUSCO short summaries in the same run have different counts",
                run_root,
                {"first": canonical, "other": other},
            )
    audit.expect(
        canonical["version"] == EXPECTED_BUSCO_VERSION,
        sample,
        stage,
        "BUSCO_VERSION_MISMATCH",
        "BUSCO version differs from the pinned workflow",
        Path(canonical["path"]),
        {"observed": canonical["version"], "expected": EXPECTED_BUSCO_VERSION},
    )
    audit.expect(
        canonical["lineage"] == EXPECTED_BUSCO_LINEAGE,
        sample,
        stage,
        "BUSCO_LINEAGE_MISMATCH",
        "BUSCO lineage is not metazoa_odb10",
        Path(canonical["path"]),
        {"observed": canonical["lineage"], "expected": EXPECTED_BUSCO_LINEAGE},
    )
    mode_ok = (
        "genome" in canonical["mode"].lower()
        if expected_mode == "genome"
        else canonical["mode"].lower() in {"protein", "proteins"}
    )
    audit.expect(
        mode_ok,
        sample,
        stage,
        "BUSCO_MODE_MISMATCH",
        f"BUSCO mode does not match expected {expected_mode} mode",
        Path(canonical["path"]),
        {"observed": canonical["mode"], "expected": expected_mode},
    )
    audit.expect(
        canonical["total"] == EXPECTED_BUSCO_GROUPS,
        sample,
        stage,
        "BUSCO_TOTAL_GROUPS_MISMATCH",
        "BUSCO total group count is not 954",
        Path(canonical["path"]),
        {"observed": canonical["total"], "expected": EXPECTED_BUSCO_GROUPS},
    )
    audit.expect(
        canonical["complete"] == canonical["single"] + canonical["duplicated"],
        sample,
        stage,
        "BUSCO_COMPLETE_ARITHMETIC_FAILURE",
        "BUSCO complete count does not equal single plus duplicated",
        Path(canonical["path"]),
        canonical,
    )
    audit.expect(
        canonical["complete"] + canonical["fragmented"] + canonical["missing"]
        == canonical["total"],
        sample,
        stage,
        "BUSCO_TOTAL_ARITHMETIC_FAILURE",
        "BUSCO complete + fragmented + missing does not equal total",
        Path(canonical["path"]),
        canonical,
    )
    for key in ("complete", "single", "duplicated", "fragmented", "missing"):
        canonical[f"{key}_percent"] = (
            100.0 * canonical[key] / canonical["total"]
            if canonical["total"]
            else 0.0
        )
    if canonical["input_path"]:
        audit.expect(
            canonical["input_path"] == str(expected_input),
            sample,
            stage,
            "BUSCO_INPUT_PATH_MISMATCH",
            "BUSCO summary input path differs from the intended input",
            Path(canonical["path"]),
            {"observed": canonical["input_path"], "expected": str(expected_input)},
        )
    else:
        audit.warning(
            sample,
            stage,
            "BUSCO_INPUT_PATH_NOT_PARSED",
            "BUSCO input path was not parsed from the short summary",
            Path(canonical["path"]),
        )

    full_table = run_root / f"run_{EXPECTED_BUSCO_LINEAGE}" / "full_table.tsv"
    if not audit.required_file(full_table, sample, stage, "BUSCO full table"):
        return dict(canonical, sample=sample, analysis=expected_mode, full_table_verified=False)
    full_table_errors_before = audit.issue_totals_by_stage[(sample, stage, "ERROR")]
    statuses_by_id: Dict[str, Set[str]] = defaultdict(set)
    sequence_ids: Set[str] = set()
    comments: List[str] = []
    with full_table.open("rt", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            if raw.startswith("#"):
                comments.append(raw.rstrip("\r\n"))
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                audit.error(
                    sample,
                    stage,
                    "BUSCO_FULL_TABLE_FIELD_COUNT",
                    f"BUSCO full table line {line_number} has fewer than two fields",
                    full_table,
                )
                continue
            busco_id, status = fields[0], fields[1]
            statuses_by_id[busco_id].add(status)
            if len(fields) >= 3 and fields[2] not in MISSING_VALUES:
                sequence_ids.add(fields[2])
    multi_status = {key: sorted(value) for key, value in statuses_by_id.items() if len(value) != 1}
    if multi_status:
        audit.error(
            sample,
            stage,
            "BUSCO_ID_MULTIPLE_STATUSES",
            "A BUSCO ID has conflicting statuses in full_table.tsv",
            full_table,
            dict(list(sorted(multi_status.items()))[:20]),
        )
    unique_status = {
        key: next(iter(value)) for key, value in statuses_by_id.items() if len(value) == 1
    }
    status_counts = Counter(unique_status.values())
    recognized = {"Complete", "Duplicated", "Fragmented", "Missing"}
    unexpected = set(status_counts) - recognized
    if unexpected:
        audit.error(
            sample,
            stage,
            "BUSCO_UNEXPECTED_STATUS",
            "Unexpected BUSCO statuses in full table",
            full_table,
            sorted(unexpected),
        )
    full_counts = {
        "single": status_counts["Complete"],
        "duplicated": status_counts["Duplicated"],
        "fragmented": status_counts["Fragmented"],
        "missing": status_counts["Missing"],
    }
    full_counts["complete"] = full_counts["single"] + full_counts["duplicated"]
    full_counts["total"] = len(unique_status)
    for key in count_keys:
        audit.expect(
            full_counts[key] == canonical[key],
            sample,
            stage,
            "BUSCO_FULL_TABLE_SUMMARY_MISMATCH",
            f"BUSCO full-table {key} count differs from short summary",
            full_table,
            {"full_table": full_counts[key], "summary": canonical[key]},
        )
    unknown_sequences = sequence_ids - allowed_sequence_ids
    if unknown_sequences:
        audit.error(
            sample,
            stage,
            "BUSCO_UNKNOWN_SEQUENCE_ID",
            "BUSCO full table contains sequence IDs absent from its declared input",
            full_table,
            {"count": len(unknown_sequences), "examples": sorted(unknown_sequences)[:20]},
        )
    comment_text = "\n".join(comments)
    lineage = re.search(r"lineage dataset is:\s*([^\s(]+)", comment_text, re.I)
    audit.expect(
        bool(lineage) and lineage.group(1) == EXPECTED_BUSCO_LINEAGE,
        sample,
        stage,
        "BUSCO_FULL_TABLE_LINEAGE_MISMATCH",
        "BUSCO full-table lineage metadata is absent or incorrect",
        full_table,
    )
    json_candidates = sorted(run_root.glob("short_summary*.json")) + sorted(
        (run_root / f"run_{EXPECTED_BUSCO_LINEAGE}").glob("short_summary*.json")
    )
    json_candidates = [
        path
        for path in json_candidates
        if path.is_file() and path.stat().st_size > 0
    ]
    if not json_candidates:
        audit.error(
            sample,
            stage,
            "MISSING_BUSCO_JSON_SUMMARY",
            "BUSCO JSON summary is missing or empty",
            run_root,
        )
    for path in json_candidates:
        audit.add_file(path, sample, stage, "BUSCO JSON summary")
        try:
            payload = json.loads(path.read_text(encoding="utf-8", errors="strict"))
        except Exception as exc:
            audit.error(
                sample,
                stage,
                "INVALID_BUSCO_JSON_SUMMARY",
                f"BUSCO JSON summary is not valid JSON: {exc}",
                path,
            )
            continue
        if not isinstance(payload, dict):
            audit.error(
                sample,
                stage,
                "INVALID_BUSCO_JSON_SUMMARY",
                "BUSCO JSON summary root is not an object",
                path,
            )
            continue
        results = payload.get("results")
        parameters = payload.get("parameters")
        lineage_data = payload.get("lineage_dataset")
        versions = payload.get("versions")
        for section_name, section in (
            ("results", results),
            ("parameters", parameters),
            ("lineage_dataset", lineage_data),
            ("versions", versions),
        ):
            audit.expect(
                isinstance(section, dict),
                sample,
                stage,
                "BUSCO_JSON_SECTION_MISSING",
                f"BUSCO JSON summary lacks an object-valued {section_name} section",
                path,
            )
        if not all(
            isinstance(section, dict)
            for section in (results, parameters, lineage_data, versions)
        ):
            continue

        json_count_keys = {
            "complete": "Complete BUSCOs",
            "single": "Single copy BUSCOs",
            "duplicated": "Multi copy BUSCOs",
            "fragmented": "Fragmented BUSCOs",
            "missing": "Missing BUSCOs",
            "total": "n_markers",
        }
        for canonical_key, json_key in json_count_keys.items():
            try:
                observed = int(results[json_key])
            except (KeyError, TypeError, ValueError):
                audit.error(
                    sample,
                    stage,
                    "BUSCO_JSON_COUNT_MISSING_OR_INVALID",
                    f"BUSCO JSON results field {json_key!r} is missing or non-integer",
                    path,
                )
                continue
            audit.expect(
                observed == canonical[canonical_key] == full_counts[canonical_key],
                sample,
                stage,
                "BUSCO_JSON_COUNT_MISMATCH",
                f"BUSCO JSON {json_key} disagrees with text/full-table counts",
                path,
                {
                    "json": observed,
                    "text": canonical[canonical_key],
                    "full_table": full_counts[canonical_key],
                },
            )

        json_lineage_name = str(lineage_data.get("name", ""))
        try:
            json_lineage_total = int(lineage_data.get("number_of_buscos"))
        except (TypeError, ValueError):
            json_lineage_total = None
        audit.expect(
            json_lineage_name == EXPECTED_BUSCO_LINEAGE
            and json_lineage_total == EXPECTED_BUSCO_GROUPS,
            sample,
            stage,
            "BUSCO_JSON_LINEAGE_MISMATCH",
            "BUSCO JSON lineage name or BUSCO-group count is incorrect",
            path,
            {
                "name": json_lineage_name,
                "number_of_buscos": json_lineage_total,
                "expected_name": EXPECTED_BUSCO_LINEAGE,
                "expected_number_of_buscos": EXPECTED_BUSCO_GROUPS,
            },
        )
        json_mode = str(parameters.get("mode", "")).lower()
        json_mode_ok = (
            "genome" in json_mode
            if expected_mode == "genome"
            else json_mode in {"protein", "proteins"}
        )
        audit.expect(
            json_mode_ok,
            sample,
            stage,
            "BUSCO_JSON_MODE_MISMATCH",
            "BUSCO JSON mode differs from the intended analysis mode",
            path,
            {"observed": json_mode, "expected": expected_mode},
        )
        json_input = str(parameters.get("in", ""))
        audit.expect(
            json_input == str(expected_input),
            sample,
            stage,
            "BUSCO_JSON_INPUT_PATH_MISMATCH",
            "BUSCO JSON input path differs from the intended input",
            path,
            {"observed": json_input, "expected": str(expected_input)},
        )
        json_lineage_parameter = Path(
            str(parameters.get("lineage_dataset", ""))
        ).name
        audit.expect(
            json_lineage_parameter == EXPECTED_BUSCO_LINEAGE,
            sample,
            stage,
            "BUSCO_JSON_LINEAGE_PARAMETER_MISMATCH",
            "BUSCO JSON lineage parameter is not metazoa_odb10",
            path,
            {
                "observed": json_lineage_parameter,
                "expected": EXPECTED_BUSCO_LINEAGE,
            },
        )
        audit.expect(
            str(versions.get("busco", "")) == EXPECTED_BUSCO_VERSION,
            sample,
            stage,
            "BUSCO_JSON_VERSION_MISMATCH",
            "BUSCO JSON version differs from the pinned workflow",
            path,
            {
                "observed": versions.get("busco"),
                "expected": EXPECTED_BUSCO_VERSION,
            },
        )
    full_table_verified = (
        audit.issue_totals_by_stage[(sample, stage, "ERROR")]
        == full_table_errors_before
    )
    return dict(
        canonical,
        sample=sample,
        analysis=expected_mode,
        full_table=str(full_table),
        full_table_verified=full_table_verified,
        full_table_sequence_ids=len(sequence_ids),
    )


def validate_swissprot_fields(
    audit: Audit,
    sample: str,
    stage: str,
    path: Path,
    line_number: int,
    fields: List[str],
    representatives: Dict[str, str],
) -> None:
    query = fields[0]
    if query not in representatives:
        audit.error(
            sample,
            stage,
            "SWISSPROT_UNKNOWN_QUERY",
            f"Swiss-Prot query {query} is absent from representative proteins",
            path,
        )
        return
    try:
        pident = float(fields[2])
        alignment_length = int(fields[3])
        query_length = int(fields[4])
        subject_length = int(fields[5])
        query_coverage = float(fields[6])
        subject_coverage = float(fields[7])
        evalue = float(fields[8])
        bitscore = float(fields[9])
    except ValueError:
        audit.error(
            sample,
            stage,
            "SWISSPROT_INVALID_NUMERIC_FIELD",
            f"Invalid numeric field on line {line_number}",
            path,
        )
        return
    expected_query_length = len(representatives[query].rstrip("*"))
    checks = [
        (0 <= pident <= 100, "pident", pident),
        (alignment_length > 0, "alignment_length", alignment_length),
        (query_length == expected_query_length, "query_length", query_length),
        (subject_length > 0, "subject_length", subject_length),
        (30 <= query_coverage <= 100, "query_coverage", query_coverage),
        (0 <= subject_coverage <= 100, "subject_coverage", subject_coverage),
        (0 <= evalue <= 1e-5, "evalue", evalue),
        (bitscore >= 0, "bitscore", bitscore),
    ]
    for passed, name, value in checks:
        if not passed:
            audit.error(
                sample,
                stage,
                "SWISSPROT_THRESHOLD_OR_LENGTH_FAILURE",
                f"Swiss-Prot {name} fails validation on line {line_number}",
                path,
                {"query": query, "value": value, "expected_query_length": expected_query_length},
            )


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "N": "ACGT",
    "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG",
    "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
}


def translate_codon(codon: str) -> str:
    if codon in CODON_TABLE:
        return CODON_TABLE[codon]
    possibilities = [IUPAC.get(base, "ACGT") for base in codon]
    amino_acids = {
        CODON_TABLE["".join(bases)] for bases in itertools.product(*possibilities)
    }
    return next(iter(amino_acids)) if len(amino_acids) == 1 else "X"


def translate_dna(sequence: str) -> str:
    return "".join(
        translate_codon(sequence[index : index + 3])
        for index in range(0, len(sequence) - 2, 3)
    )


def flatten_metrics(sample: str, metrics: Dict[str, object]) -> dict:
    row = {"sample": sample, "label": SAMPLES[sample]["label"]}
    row.update(metrics)
    return row


def stage_key_metrics(metrics: Dict[str, object], stage: str) -> Dict[str, object]:
    """Select compact, stage-relevant metrics for the PASS/WARN/FAIL ledger."""
    prefixes = {
        "input": ("input_", "locked_source_", "short_fasta_"),
        "repeatmasker": ("repeatmasker_",),
        "genome_busco": ("genome_busco_",),
        "braker": ("braker_",),
        "protein_busco": (
            "protein_busco_",
            "all_protein_",
            "representative_",
            "represented_genes",
        ),
        "eggnog": ("eggnog_",),
        "swissprot": ("swissprot_",),
        "interpro": ("interpro",),
    }[stage]
    return {
        key: metrics[key]
        for key in sorted(metrics)
        if any(key.startswith(prefix) for prefix in prefixes)
    }


def build_stage_rows(audit: Audit) -> List[dict]:
    rows: List[dict] = []
    for sample in SAMPLES:
        for stage in STAGES:
            counts = {
                severity: audit.issue_totals_by_stage[(sample, stage, severity)]
                for severity in ("ERROR", "WARNING", "INFO")
            }
            status = (
                "FAIL"
                if counts["ERROR"]
                else "WARN"
                if counts["WARNING"]
                else "PASS"
            )
            file_count = sum(
                1
                for record in audit.files.values()
                if record.sample == sample and record.stage == stage
            )
            stage_files = sorted(
                (
                    {
                        "role": record.role,
                        "path": record.path,
                        "sha256": record.sha256,
                    }
                    for record in audit.files.values()
                    if record.sample == sample and record.stage == stage
                ),
                key=lambda row: (row["path"], row["role"]),
            )
            rows.append(
                {
                    "sample": sample,
                    "label": SAMPLES[sample]["label"],
                    "stage": stage,
                    "status": status,
                    "attempted": (sample, stage) in audit.attempted_stages,
                    "validated_without_error": (
                        (sample, stage) in audit.attempted_stages
                        and counts["ERROR"] == 0
                    ),
                    "error_count": counts["ERROR"],
                    "warning_count": counts["WARNING"],
                    "info_count": counts["INFO"],
                    "audited_file_count": file_count,
                    "key_metrics": stage_key_metrics(audit.metrics[sample], stage),
                    "file_sha256": stage_files,
                }
            )
        sample_rows = [row for row in rows if row["sample"] == sample]
        rows.append(
            {
                "sample": sample,
                "label": SAMPLES[sample]["label"],
                "stage": "ALL",
                "status": (
                    "FAIL"
                    if any(row["status"] == "FAIL" for row in sample_rows)
                    else "WARN"
                    if any(row["status"] == "WARN" for row in sample_rows)
                    else "PASS"
                ),
                "attempted": all(row["attempted"] for row in sample_rows),
                "validated_without_error": all(
                    row["validated_without_error"] for row in sample_rows
                ),
                "error_count": sum(row["error_count"] for row in sample_rows),
                "warning_count": sum(row["warning_count"] for row in sample_rows),
                "info_count": sum(row["info_count"] for row in sample_rows),
                "audited_file_count": sum(
                    row["audited_file_count"] for row in sample_rows
                ),
                "key_metrics": {},
                "file_sha256": [],
            }
        )
    sample_overall = [row for row in rows if row["stage"] == "ALL"]
    rows.append(
        {
            "sample": "ALL",
            "label": "LHT-04 and NJ-lan",
            "stage": "ALL",
            "status": (
                "FAIL"
                if any(row["status"] == "FAIL" for row in sample_overall)
                else "WARN"
                if any(row["status"] == "WARN" for row in sample_overall)
                else "PASS"
            ),
            "attempted": all(row["attempted"] for row in sample_overall),
            "validated_without_error": all(
                row["validated_without_error"] for row in sample_overall
            ),
            "error_count": sum(row["error_count"] for row in sample_overall),
            "warning_count": sum(row["warning_count"] for row in sample_overall),
            "info_count": sum(row["info_count"] for row in sample_overall),
            "audited_file_count": len(audit.files),
            "key_metrics": {},
            "file_sha256": [],
        }
    )
    return rows


def write_reports(audit: Audit, output_dir: Path, started_at: str) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    status = (
        "FAIL"
        if audit.issue_totals["ERROR"]
        else "WARN"
        if audit.issue_totals["WARNING"]
        else "PASS"
    )
    finished_at = datetime.now(timezone.utc).isoformat()
    file_rows = [asdict(record) for record in sorted(audit.files.values(), key=lambda x: x.path)]
    stage_rows = build_stage_rows(audit)
    interpretation_notes = [
        "A PASS means the fixed files and cross-file relationships tested here were internally consistent; it is not biological validation of every prediction.",
        "BRAKER used compleasm/BUSCO lineage hints and compleasm-based model selection, so protein BUSCO is not fully independent of gene prediction.",
        "samples.csv declares protein evidence but no RNA-seq or Iso-Seq evidence; RNA-labelled support columns must not be interpreted as experimental RNA support.",
        "eggNOG and InterProScan may legitimately omit representative proteins with no hit; their query IDs must nevertheless be subsets of the representative FASTA.",
    ]
    report = {
        "schema_version": 1,
        "audit_script_version": AUDIT_VERSION,
        "status": status,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "host": socket.gethostname(),
        "root": str(audit.root),
        "python": sys.version,
        "expected": {
            "samples": SAMPLES,
            "busco_version": EXPECTED_BUSCO_VERSION,
            "busco_lineage": EXPECTED_BUSCO_LINEAGE,
            "busco_groups": EXPECTED_BUSCO_GROUPS,
        },
        "metrics": audit.metrics,
        "busco": audit.busco_rows,
        "functional_annotation": audit.functional_rows,
        "stage_results": stage_rows,
        "issue_counts": dict(audit.issue_totals),
        "issue_counts_by_stage": [
            {
                "sample": sample,
                "stage": stage,
                "severity": severity,
                "count": count,
            }
            for (sample, stage, severity), count in sorted(
                audit.issue_totals_by_stage.items()
            )
        ],
        "suppressed_issue_counts": dict(audit.suppressed_issue_totals),
        "issues": [asdict(issue) for issue in audit.issues],
        "audited_files": file_rows,
        "interpretation_notes": interpretation_notes,
    }
    json_path = output_dir / "final_annotation_audit.json"
    atomic_write_text(json_path, json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    metric_rows = [flatten_metrics(sample, audit.metrics[sample]) for sample in SAMPLES]
    metric_columns = ["sample", "label"] + sorted(
        {key for row in metric_rows for key in row if key not in {"sample", "label"}}
    )
    atomic_write_text(
        output_dir / "final_annotation_stats.tsv", tsv_text(metric_rows, metric_columns)
    )
    busco_columns = [
        "sample", "analysis", "version", "lineage", "mode", "input_path",
        "complete", "single", "duplicated", "fragmented", "missing", "total",
        "complete_percent", "single_percent", "duplicated_percent",
        "fragmented_percent", "missing_percent",
        "full_table_verified", "path", "full_table",
    ]
    atomic_write_text(
        output_dir / "busco_summary.tsv", tsv_text(audit.busco_rows, busco_columns)
    )
    functional_columns = [
        "sample", "source", "data_rows", "unique_queries",
        "representative_proteins", "annotated_queries", "annotated_percent",
        "details",
    ]
    atomic_write_text(
        output_dir / "functional_annotation_summary.tsv",
        tsv_text(audit.functional_rows, functional_columns),
    )
    issue_rows = [asdict(issue) for issue in audit.issues]
    issue_columns = ["severity", "sample", "stage", "code", "message", "path", "details"]
    atomic_write_text(
        output_dir / "final_annotation_issues.tsv", tsv_text(issue_rows, issue_columns)
    )
    file_columns = [
        "sample", "stage", "role", "path", "bytes_before", "mtime_ns_before",
        "sha256", "stable_during_audit",
    ]
    atomic_write_text(
        output_dir / "audited_file_manifest.tsv", tsv_text(file_rows, file_columns)
    )

    stage_columns = [
        "sample",
        "label",
        "stage",
        "status",
        "attempted",
        "validated_without_error",
        "error_count",
        "warning_count",
        "info_count",
        "audited_file_count",
        "key_metrics",
        "file_sha256",
    ]
    stage_tsv = tsv_text(stage_rows, stage_columns)
    atomic_write_text(output_dir / "results_audit.tsv", stage_tsv)

    final_summary = {
        "schema_version": 1,
        "status": status,
        "generated_at_utc": finished_at,
        "root": str(audit.root),
        "samples": audit.metrics,
        "busco": audit.busco_rows,
        "functional_annotation": audit.functional_rows,
        "detailed_audit": str(json_path),
    }
    results_audit = {
        "schema_version": 1,
        "status": status,
        "generated_at_utc": finished_at,
        "root": str(audit.root),
        "stage_results": stage_rows,
        "metrics": audit.metrics,
        "issue_counts": dict(audit.issue_totals),
        "suppressed_issue_counts": dict(audit.suppressed_issue_totals),
        "issues": [asdict(issue) for issue in audit.issues],
        "audited_files": file_rows,
        "detailed_audit": str(json_path),
        "interpretation_notes": interpretation_notes,
    }
    final_summary_json = json.dumps(
        final_summary, indent=2, ensure_ascii=False
    ) + "\n"
    results_audit_json = json.dumps(
        results_audit, indent=2, ensure_ascii=False
    ) + "\n"
    final_summary_tsv = tsv_text(metric_rows, metric_columns)

    # Detailed copies live under ROOT/08_summary; standardized package gates are
    # the four top-level files required by PACKAGE_PLAN.md.
    atomic_write_text(output_dir / "final_summary.tsv", final_summary_tsv)
    atomic_write_text(output_dir / "final_summary.json", final_summary_json)
    atomic_write_text(output_dir / "results_audit.json", results_audit_json)
    atomic_write_text(audit.root / "final_summary.tsv", final_summary_tsv)
    atomic_write_text(audit.root / "final_summary.json", final_summary_json)
    atomic_write_text(audit.root / "results_audit.tsv", stage_tsv)
    atomic_write_text(audit.root / "results_audit.json", results_audit_json)

    # Each sample gets a compact self-contained view in its pre-created
    # 08_summary directory.  These are useful when a sample is reviewed alone.
    for sample in SAMPLES:
        sample_dir = audit.root / sample / "08_summary"
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_metric_rows = [
            row for row in metric_rows if row["sample"] == sample
        ]
        sample_busco = [row for row in audit.busco_rows if row["sample"] == sample]
        sample_functional = [
            row for row in audit.functional_rows if row["sample"] == sample
        ]
        sample_stage_rows = [
            row for row in stage_rows if row["sample"] == sample
        ]
        sample_issues = [
            asdict(issue) for issue in audit.issues if issue.sample == sample
        ]
        sample_files = [row for row in file_rows if row["sample"] == sample]
        sample_status = (
            "FAIL"
            if any(row["status"] == "FAIL" for row in sample_stage_rows)
            else "WARN"
            if any(row["status"] == "WARN" for row in sample_stage_rows)
            else "PASS"
        )
        sample_summary = {
            "schema_version": 1,
            "sample": sample,
            "label": SAMPLES[sample]["label"],
            "status": sample_status,
            "generated_at_utc": finished_at,
            "metrics": audit.metrics[sample],
            "busco": sample_busco,
            "functional_annotation": sample_functional,
        }
        sample_audit = {
            "schema_version": 1,
            "sample": sample,
            "label": SAMPLES[sample]["label"],
            "status": sample_status,
            "generated_at_utc": finished_at,
            "stage_results": sample_stage_rows,
            "issues": sample_issues,
            "audited_files": sample_files,
            "interpretation_notes": interpretation_notes,
        }
        atomic_write_text(
            sample_dir / "final_summary.tsv",
            tsv_text(sample_metric_rows, metric_columns),
        )
        atomic_write_text(
            sample_dir / "final_summary.json",
            json.dumps(sample_summary, indent=2, ensure_ascii=False) + "\n",
        )
        atomic_write_text(
            sample_dir / "results_audit.tsv",
            tsv_text(sample_stage_rows, stage_columns),
        )
        atomic_write_text(
            sample_dir / "results_audit.json",
            json.dumps(sample_audit, indent=2, ensure_ascii=False) + "\n",
        )
        atomic_write_text(
            sample_dir / "busco_summary.tsv",
            tsv_text(sample_busco, busco_columns),
        )
        atomic_write_text(
            sample_dir / "functional_annotation_summary.tsv",
            tsv_text(sample_functional, functional_columns),
        )
        atomic_write_text(
            sample_dir / "audited_file_manifest.tsv",
            tsv_text(sample_files, file_columns),
        )
    return report


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        help="Pipeline root, e.g. /work1/sponge/lht-njlan-final-annotation-20261007",
    )
    parser.add_argument(
        "--max-issues",
        type=int,
        default=2000,
        help="Maximum detailed issues retained in JSON/TSV (default: 2000)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        print(f"ERROR: pipeline root is not a directory: {root}", file=sys.stderr)
        return 2
    if args.max_issues < 1:
        print("ERROR: --max-issues must be positive", file=sys.stderr)
        return 2
    started_at = datetime.now(timezone.utc).isoformat()
    audit = Audit(root, max_issues=args.max_issues)
    for sample, config in SAMPLES.items():
        SampleAuditor(audit, sample, config).run()
    audit.finalize_file_hashes()
    report = write_reports(audit, root / "08_summary", started_at)
    print(
        json.dumps(
            {
                "status": report["status"],
                "errors": audit.issue_totals["ERROR"],
                "warnings": audit.issue_totals["WARNING"],
                "report": str(root / "08_summary" / "final_annotation_audit.json"),
            },
            ensure_ascii=False,
        )
    )
    return 1 if report["status"] == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
