#!/usr/bin/env python3
"""Add explicit gene= tags to audited longest-isoform FASTA headers.

The source FASTAs already contain one representative transcript per gene and
use BRAKER identifiers of the form g<number>.t<number>.  This script preserves
the sequences and protein identifiers verbatim while adding the gene identifier
required by the reusable functional-annotation integration workflow.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


PROTEIN_ID = re.compile(r"^(g\d+)\.t\d+$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize(source: Path, destination: Path) -> dict:
    records = 0
    genes: set[str] = set()
    sequence_lines = 0
    with source.open(encoding="utf-8") as src, destination.open("w", encoding="utf-8") as dst:
        for line_number, raw in enumerate(src, 1):
            if raw.startswith(">"):
                header = raw[1:].rstrip("\r\n")
                protein_id = header.split()[0]
                match = PROTEIN_ID.fullmatch(protein_id)
                if match is None:
                    raise ValueError(f"{source}:{line_number}: unexpected protein identifier {protein_id!r}")
                gene_id = match.group(1)
                if gene_id in genes:
                    raise ValueError(f"{source}:{line_number}: duplicate representative for {gene_id}")
                genes.add(gene_id)
                records += 1
                dst.write(f">{header} gene={gene_id}\n")
            else:
                if not raw.strip():
                    continue
                sequence_lines += 1
                dst.write(raw if raw.endswith("\n") else raw + "\n")
    if records == 0 or sequence_lines == 0:
        raise ValueError(f"{source}: empty FASTA")
    return {
        "source": str(source),
        "destination": str(destination),
        "records": records,
        "unique_genes": len(genes),
        "source_sha256": sha256(source),
        "destination_sha256": sha256(destination),
        "transformation": "append gene=<BRAKER gene ID> derived from gN.tM; sequence and protein ID unchanged",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", nargs="+", required=True, metavar="SOURCE=DESTINATION")
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()
    audit = []
    for pair in args.pairs:
        source_text, destination_text = pair.split("=", 1)
        source = Path(source_text).resolve()
        destination = Path(destination_text).resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        audit.append(normalize(source, destination))
    Path(args.audit).write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
