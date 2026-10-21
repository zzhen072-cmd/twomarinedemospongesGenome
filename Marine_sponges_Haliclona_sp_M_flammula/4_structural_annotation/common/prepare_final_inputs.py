#!/usr/bin/env python3
"""Create audited short-ID FASTA inputs for final LHT-04 and NJ-lan assemblies."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterator, TextIO, Tuple


SAMPLES = {
    "lht04": {
        "label": "LHT-04",
        "prefix": "lhtf_ctg",
        "source": "/work1/sponge/lht_reassembly_compare/15_read_validated_curated/LHT_curated_primary.fasta",
        "expected_records": 844,
        "expected_bp": 115_734_075,
        "expected_sha256": "026d8b2e3e9888de0763d42ce226fcf7c63062a145e6ad8478e2e8c1725b2499",
        "min_contig": 10_000,
    },
    "njlan": {
        "label": "NJ-lan",
        "prefix": "njlf_ctg",
        "source": "/work1/sponge/lan_tier2_decontam_review/assembly/lan_host_review_filtered.fasta",
        "expected_records": 7_117,
        "expected_bp": 291_889_606,
        "expected_sha256": "cea9b34517164f02103111a3a1dfa7cfd65ec95a61c224da1901036a4cca120a",
        "min_contig": 2_000,
    },
}

PROTEIN_EVIDENCE = "/work1/sponge/genome_annotation/protein_evidence/Porifera_2species.faa"
BUSCO_DOWNLOADS = "/work1/sponge/genome_annotation/braker4_shared/busco_downloads"
COMPLEASM_DOWNLOADS = "/work1/sponge/genome_annotation/braker4_shared/compleasm_downloads"
BRAKER_IMAGE = "/work1/sponge/genome_annotation/braker4_containers/braker3_v3.1.1.sif"
AGAT_IMAGE = "/work1/sponge/genome_annotation/braker4_containers/agat_1.4.1.sif"


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("rt", encoding="ascii")


def open_binary(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def uncompressed_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open_binary(path) as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    name = None
    chunks = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = raw[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name = header.split()[0]
                chunks = []
            else:
                if name is None:
                    if raw.strip():
                        raise ValueError(f"{path}:{line_number}: sequence before header")
                    continue
                chunks.append("".join(raw.split()))
    if name is not None:
        yield name, "".join(chunks)


def write_wrapped(handle: TextIO, sequence: str, width: int = 60) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def braker_config(root: Path, sample: str, min_contig: int) -> str:
    run = root / sample / "03_braker4"
    return f"""[paths]
samples_file = {run}/samples.csv
augustus_config_path = {run}/augustus_config
busco_download_path = {BUSCO_DOWNLOADS}
compleasm_download_path = {COMPLEASM_DOWNLOADS}

[containers]
braker3_image = {BRAKER_IMAGE}
agat_image = {AGAT_IMAGE}

[PARAMS]
fungus = 0
min_contig = {min_contig}
use_varus = 0
skip_optimize_augustus = 1
skip_single_exon_downsampling = 0
downsampling_lambda = 2
use_dev_shm = 0
use_compleasm_hints = 1
skip_busco = 1
run_omark = 0
translation_table = 1
gc_donor = 0.001
allow_hinted_splicesites = gcag,atac
augustus_chunksize = 3000000
augustus_overlap = 500000
run_ncrna = 0
run_best_by_compleasm = 1
masking_tool = repeatmasker
use_minisplice = 0
no_cleanup = 0

[fantasia]
enable = 0

[SLURM_ARGS]
cpus_per_task = 30
mem_of_node = 108000
max_runtime = 1440
"""


def samples_csv(root: Path, sample: str) -> str:
    sample_name = f"{sample}_final"
    fasta = root / sample / "00_input" / f"{sample_name}.shortids.fasta"
    masked = root / sample / "02_repeatmasker" / f"{sample_name}.shortids.fasta.masked"
    header = (
        "sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,"
        "fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,"
        "busco_lineage,reference_gtf\n"
    )
    row = f"{sample_name},{fasta},{masked},{PROTEIN_EVIDENCE},,,,,,,,,metazoa_odb10,\n"
    return header + row


def prepare_one(root: Path, sample: str, config: dict) -> dict:
    source = Path(config["source"])
    if not source.is_file():
        gz = Path(str(source) + ".gz")
        if not gz.is_file():
            raise FileNotFoundError(source)
        source = gz
    observed_sha = uncompressed_sha256(source)
    if observed_sha != config["expected_sha256"]:
        raise ValueError(
            f"{sample}: input SHA-256 {observed_sha} != {config['expected_sha256']}"
        )

    base = root / sample
    for directory in (
        "00_input",
        "01_genome_busco",
        "02_repeatmasker",
        "03_braker4",
        "04_protein_busco",
        "05_eggnog",
        "06_swissprot",
        "07_interpro",
        "08_summary",
        "logs",
    ):
        (base / directory).mkdir(parents=True, exist_ok=True)

    sample_name = f"{sample}_final"
    fasta_out = base / "00_input" / f"{sample_name}.shortids.fasta"
    map_out = base / "00_input" / "contig_id_map.tsv"
    manifest_out = base / "00_input" / "input_manifest.json"
    for target in (fasta_out, map_out, manifest_out):
        if target.exists():
            raise FileExistsError(f"refusing to overwrite {target}")

    seen = set()
    record_count = 0
    total_bp = 0
    sequence_digest = hashlib.sha256()
    with fasta_out.open("wt", encoding="ascii") as fasta_handle, map_out.open(
        "wt", encoding="utf-8"
    ) as map_handle:
        map_handle.write("short_id\toriginal_id\tlength\tsequence_sha256\n")
        for record_count, (original_id, sequence) in enumerate(fasta_records(source), 1):
            if original_id in seen:
                raise ValueError(f"{sample}: duplicate FASTA ID {original_id}")
            seen.add(original_id)
            if not sequence:
                raise ValueError(f"{sample}: empty sequence {original_id}")
            if set(sequence.upper()) - set("ACGTNRYKMSWBDHV"):
                raise ValueError(f"{sample}: unsupported sequence character in {original_id}")
            short_id = f"{config['prefix']}{record_count:06d}"
            sequence_hash = hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()
            sequence_digest.update(short_id.encode("ascii"))
            sequence_digest.update(b"\0")
            sequence_digest.update(sequence.upper().encode("ascii"))
            sequence_digest.update(b"\0")
            total_bp += len(sequence)
            fasta_handle.write(f">{short_id}\n")
            write_wrapped(fasta_handle, sequence)
            map_handle.write(
                f"{short_id}\t{original_id}\t{len(sequence)}\t{sequence_hash}\n"
            )

    if record_count != config["expected_records"] or total_bp != config["expected_bp"]:
        raise ValueError(
            f"{sample}: observed {record_count} records/{total_bp} bp; expected "
            f"{config['expected_records']} records/{config['expected_bp']} bp"
        )

    manifest = {
        "sample": config["label"],
        "source_path": str(source),
        "source_uncompressed_sha256": observed_sha,
        "short_fasta_path": str(fasta_out),
        "short_fasta_sha256": hashlib.sha256(fasta_out.read_bytes()).hexdigest(),
        "contig_id_map": str(map_out),
        "records": record_count,
        "bp": total_bp,
        "ordered_shortid_sequence_sha256": sequence_digest.hexdigest(),
        "status": "PASS",
    }
    manifest_out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    run = base / "03_braker4"
    (run / "config.ini").write_text(
        braker_config(root, sample, int(config["min_contig"])), encoding="utf-8"
    )
    (run / "samples.csv").write_text(samples_csv(root, sample), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", default="/work1/sponge/lht-njlan-final-annotation-20261007"
    )
    args = parser.parse_args()
    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    manifests = {sample: prepare_one(root, sample, cfg) for sample, cfg in SAMPLES.items()}
    combined = root / "input_audit.json"
    if combined.exists():
        raise FileExistsError(f"refusing to overwrite {combined}")
    combined.write_text(
        json.dumps({"status": "PASS", "samples": manifests}, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / ".inputs_prepared").write_text("PASS\n", encoding="ascii")
    print(json.dumps({"status": "PASS", "root": str(root), "samples": manifests}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
