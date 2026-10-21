#!/usr/bin/env python3
"""Reproducible, gene-level functional annotation tables (Python 3.9+; stdlib).

Manifest example (paths are relative to the manifest, or absolute):
{
  "species": {
    "gzx": {"proteins": "gzx.aa", "eggnog": "gzx.annotations",
            "interpro": "gzx.interproscan.tsv", "swissprot": "gzx.best.tsv"},
    "S04": {"proteins": "S04.aa", "eggnog": "S04.annotations",
            "interpro": "S04.interproscan.tsv", "swissprot": null}
  },
  "go": "go-basic.obo",
  "kegg": "ko00001.keg"
}

Every FASTA header must contain gene=, with exactly one protein per unique gene.
No ID conversion, significance threshold, or best-hit selection is invented.
SwissProt without a declared schema uses only the first field as query ID and
preserves the original row. It must already be the desired filtered/best-hit file.

Any input may instead be {"path": "..."}. SwissProt additionally supports
"columns": ["qseqid", ...], "has_header": false, "query_column": "qseqid",
and "functional_columns": {"go": "GO_column", "ko": "KO_column",
"pathway": "pathway_column", "cog": "COG_column"}. Functional columns are
used only when explicitly declared. A real header can be read with has_header.

Default UpSet sources: eggNOG seed_ortholog hit, real InterPro accession,
SwissProt input hit. Missing source files are NA, never zero. The InterProScan
any-signature and Pfam matches are independent columns/coverage metrics.

Primary KEGG counts use the union of directly supplied pathway annotations.
KO-derived pathway IDs are exported separately and do not inflate those counts.
All percentages use the number of representative FASTA genes as denominator.
GO bins are all matching is_a/part_of generic-slim ancestors, excluding roots.

On validation failure, only audit reports and analysis_manifest.json are written;
functional results are not generated. Read manifest status before using tables.
An output directory containing existing functional result tables is refused to
prevent stale results from being mistaken for a successful new analysis.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import itertools
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from ontology import GOTerm, GOOntology, load_go, load_kegg, normalize_pathway, sha256_file

GO_RE = re.compile(r"(?<![A-Za-z0-9_])GO:\d{7}(?!\d)")
KO_RE = re.compile(r"(?<![A-Za-z0-9_])K\d{5}(?!\d)")
IPR_RE = re.compile(r"(?<![A-Za-z0-9_])IPR\d{6}(?!\d)")
PFAM_RE = re.compile(r"(?<![A-Za-z0-9_])PF\d{5}(?:\.\d+)?(?!\d)")
GENE_RE = re.compile(r"(?:^|[\s;\[])gene=([^\s;\]]+)")
MISSING = {"", "-", ".", "NA", "N/A", "null", "None"}
SOURCE_ORDER = ("eggNOG", "InterPro", "SwissProt")
SOURCE_KEYS = {"eggNOG": "eggnog", "InterPro": "interpro", "SwissProt": "swissprot"}
ONTOLOGY_LABELS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
IPR_COLUMNS = (
    "protein_accession", "sequence_md5", "sequence_length", "analysis",
    "signature_accession", "signature_description", "start", "stop", "score",
    "status", "date", "interpro_accession", "interpro_description", "GO_annotations", "pathways",
)
EGG_REQUIRED = ("query", "seed_ortholog", "COG_category", "GOs", "KEGG_ko", "KEGG_Pathway")
EGG_EXPORT = (
    "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "Description",
    "Preferred_name", "EC", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
)
COG_NAMES = {
    "A": "RNA processing and modification", "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion", "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism", "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only", "S": "Function unknown", "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport", "V": "Defense mechanisms",
    "W": "Extracellular structures", "X": "Mobilome: prophages, transposons", "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}
FUNCTIONAL_RESULTS = (
    "gene_functional_annotations.tsv", "functional_coverage.tsv", "upset_intersections.tsv",
    "upset_set_sizes.tsv", "goslim_counts.tsv", "kegg_level2_counts.tsv",
    "kegg_pathway_counts.tsv", "protein_lengths.tsv", "protein_length_summary.tsv", "cog_counts.tsv",
)


def is_present(value):
    return str(value).strip() not in MISSING


def json_text(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def joined(values):
    return ";".join(sorted(set(values)))


def percent(count, total):
    return round(100.0 * count / total, 6) if total else 0.0


def flag(value, available=True):
    return ("true" if value else "false") if available else "NA"


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else Path(path).open(encoding="utf-8")


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class Audit:
    def __init__(self):
        self.issues = []
        self.query_counts = defaultdict(Counter)
        self.query_first_line = {}
        self.source_info = {}
        self.length_checks = Counter()

    def issue(self, species, source, kind, detail, line="", query_id="", severity="error"):
        self.issues.append(dict(species=species, source=source, severity=severity,
                                kind=kind, line=line, query_id=query_id, detail=str(detail)))

    def query(self, species, source, query_id, line):
        self.query_counts[(species, source)][query_id] += 1
        self.query_first_line.setdefault((species, source, query_id), line)

    def errors(self):
        return [item for item in self.issues if item["severity"] == "error"]


def read_fasta(handle, species, audit):
    proteins = {}
    genes = {}
    current = None
    sequence = []

    def flush():
        if current is None:
            return
        query_id, gene_id, header, start_line = current
        seq = "".join(sequence).upper()
        if not seq:
            audit.issue(species, "proteins", "empty_sequence", "Protein sequence is empty", start_line, query_id)
        invalid = sorted(set(seq) - set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"))
        if invalid:
            audit.issue(species, "proteins", "invalid_sequence_symbols", joined(invalid), start_line, query_id)
        terminal_stops = len(seq) - len(seq.rstrip("*"))
        internal_stops = seq.count("*") - terminal_stops
        if internal_stops:
            audit.issue(species, "proteins", "internal_stop", f"{internal_stops} internal stop markers", start_line, query_id, "warning")
        if query_id in proteins:
            audit.issue(species, "proteins", "duplicate_protein_id", "FASTA protein IDs must be unique", start_line, query_id)
            return
        if gene_id in genes:
            audit.issue(species, "proteins", "multiple_representatives_per_gene", f"gene={gene_id}; previous protein={genes[gene_id]}", start_line, query_id)
        genes[gene_id] = query_id
        proteins[query_id] = {
            "gene_id": gene_id, "protein_id": query_id, "fasta_header": header,
            "protein_aa": len(seq.replace("*", "")), "sequence_characters": len(seq),
            "terminal_stop_count": terminal_stops, "internal_stop_count": internal_stops,
        }

    for line_number, raw in enumerate(handle, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:]
            pieces = header.split()
            query_id = pieces[0] if pieces else ""
            gene_ids = {value.strip("\"'") for value in GENE_RE.findall(header)}
            if not query_id or len(gene_ids) != 1:
                audit.issue(species, "proteins", "invalid_fasta_header", "Require a protein ID and exactly one distinct gene= value", line_number, query_id)
            gene_id = next(iter(gene_ids)) if len(gene_ids) == 1 else ""
            current = (query_id, gene_id, header, line_number)
            sequence = []
        elif current is None:
            audit.issue(species, "proteins", "sequence_before_header", "Sequence encountered before first FASTA header", line_number)
        else:
            sequence.append("".join(line.split()))
    flush()
    if not proteins:
        audit.issue(species, "proteins", "empty_fasta", "No protein records found")
    return proteins


def read_eggnog(handle, species, proteins, audit):
    records = defaultdict(list)
    header = None
    comments = []
    for line_number, raw in enumerate(handle, 1):
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("#query\t") or line == "#query":
            current_header = line[1:].split("\t")
            if header is not None and header != current_header:
                audit.issue(species, "eggNOG", "inconsistent_header", "Multiple different #query headers", line_number)
            header = current_header
            if len(header) != len(set(header)):
                audit.issue(species, "eggNOG", "duplicate_column_name", json_text(header), line_number)
            absent = set(EGG_REQUIRED) - set(header)
            if absent:
                audit.issue(species, "eggNOG", "missing_required_columns", joined(absent), line_number)
            continue
        if line.startswith("#"):
            comments.append(line)
            continue
        fields = line.split("\t")
        query_id = fields[0]
        audit.query(species, "eggNOG", query_id, line_number)
        if header is None:
            audit.issue(species, "eggNOG", "data_before_query_header", "A real #query header is required; positions are never assumed", line_number, query_id)
            continue
        if len(fields) != len(header):
            audit.issue(species, "eggNOG", "row_width_mismatch", f"Expected {len(header)} fields, found {len(fields)}", line_number, query_id)
            continue
        row = dict(zip(header, fields))
        query_id = row["query"]
        if query_id in proteins:
            records[query_id].append(row)
    if header is None:
        audit.issue(species, "eggNOG", "missing_query_header", "No real #query header found")
    audit.source_info[(species, "eggNOG")] = {"header": header, "comments": comments}
    return records


def read_interpro(handle, species, proteins, audit):
    # Reactome annotations can make this TSV several GB. Keep compact domains,
    # term sets and original line references, never the huge raw pathway strings.
    records = {}
    widths = Counter()
    mismatch_seen = set()
    for line_number, raw in enumerate(handle, 1):
        line = raw.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        query_id = fields[0]
        audit.query(species, "InterPro", query_id, line_number)
        widths[len(fields)] += 1
        if len(fields) < 11:
            audit.issue(species, "InterPro", "short_interpro_row", "Standard InterProScan TSV requires at least 11 fields", line_number, query_id)
            continue
        row = {IPR_COLUMNS[index] if index < len(IPR_COLUMNS) else f"extra_{index+1}": value
               for index, value in enumerate(fields)}
        try:
            reported_length = int(fields[2])
            start, stop = int(fields[6]), int(fields[7])
            if not 1 <= start <= stop <= reported_length:
                raise ValueError("domain coordinates outside sequence")
        except ValueError as exc:
            audit.issue(species, "InterPro", "invalid_length_or_coordinates", str(exc), line_number, query_id)
            continue
        if query_id not in proteins:
            continue
        protein = proteins[query_id]
        exact = reported_length == protein["protein_aa"]
        terminal_stop_tolerance = (
            protein["terminal_stop_count"] == 1 and protein["internal_stop_count"] == 0
            and reported_length == protein["protein_aa"] + 1
        )
        comparison = "exact" if exact else "terminal_stop_plus_one" if terminal_stop_tolerance else "mismatch"
        audit.length_checks[(species, comparison)] += 1
        row["_fasta_length_comparison"] = comparison
        if not exact and not terminal_stop_tolerance and (query_id, reported_length) not in mismatch_seen:
            audit.issue(species, "InterPro", "protein_length_mismatch",
                        f"InterProScan length={reported_length}; FASTA aa={protein['protein_aa']}; terminal_stops={protein['terminal_stop_count']}; only a known single terminal '*' permits +1",
                        line_number, query_id)
            mismatch_seen.add((query_id, reported_length))
        if query_id not in records:
            records[query_id] = dict(row_count=0, line_numbers=[], domains=[], go=set(), ko=set(),
                                    pathway_tokens=set(), cog=set(), ipr_ids=set(), signatures=set(),
                                    pfam_ids=set(), analyses=set())
        aggregate = records[query_id]
        aggregate["row_count"] += 1
        aggregate["line_numbers"].append(line_number)
        # Columns 2-14 retain sequence checksum, analysis, exact signature,
        # coordinates, score, status/date, InterPro and GO evidence text.
        aggregate["domains"].append(fields[1:14])
        aggregate["go"].update(GO_RE.findall(row.get("GO_annotations", "")))
        aggregate["pathway_tokens"].update(pathway_tokens(row.get("pathways", ""), interpro=True))
        aggregate["ipr_ids"].update(IPR_RE.findall(row.get("interpro_accession", "")))
        aggregate["analyses"].add(row["analysis"])
        if is_present(row["signature_accession"]):
            aggregate["signatures"].add(f"{row['analysis']}:{row['signature_accession']}")
            if row["analysis"].casefold() == "pfam":
                aggregate["pfam_ids"].add(row["signature_accession"])
    audit.source_info[(species, "InterPro")] = {
        "standard_columns": list(IPR_COLUMNS), "row_width_counts": dict(widths),
        "compact_domain_columns": list(IPR_COLUMNS[1:14]),
        "raw_storage": "Original input file retained; table provides one-based raw line numbers. Giant pathway strings are not copied into memory or the gene table.",
    }
    return records


def read_swissprot(handle, species, proteins, audit, spec):
    records = defaultdict(list)
    columns = spec.get("columns")
    if columns is not None and (not isinstance(columns, list) or not all(isinstance(value, str) for value in columns)):
        raise ValueError("SwissProt columns must be a list of field names")
    has_header = bool(spec.get("has_header", False))
    consumed_header = False
    comments = []
    functional_columns = spec.get("functional_columns", {})
    for line_number, raw in enumerate(handle, 1):
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("#"):
            comments.append(line)
            continue
        fields = line.split("\t")
        if has_header and not consumed_header:
            consumed_header = True
            if columns is not None and columns != fields:
                audit.issue(species, "SwissProt", "declared_header_mismatch", "Actual header differs from declared columns", line_number)
            columns = fields
            continue
        if columns is not None:
            if len(columns) != len(set(columns)):
                audit.issue(species, "SwissProt", "duplicate_column_name", json_text(columns), line_number)
                continue
            query_column = spec.get("query_column", columns[0])
            if query_column not in columns:
                raise ValueError(f"SwissProt query_column {query_column!r} is absent from columns")
            query_position = columns.index(query_column)
            query_id = fields[query_position] if query_position < len(fields) else fields[0]
            audit.query(species, "SwissProt", query_id, line_number)
            if len(fields) != len(columns):
                audit.issue(species, "SwissProt", "row_width_mismatch", f"Expected {len(columns)}, found {len(fields)}", line_number, query_id)
                continue
            row = dict(zip(columns, fields))
        else:
            query_id = fields[0]
            audit.query(species, "SwissProt", query_id, line_number)
            row = {"query": query_id, "_raw_fields": fields}
        row["_raw_tsv"] = line
        if query_id in proteins:
            records[query_id].append(row)
    if has_header and not consumed_header:
        audit.issue(species, "SwissProt", "missing_declared_header", "has_header=true but no header was found")
    if functional_columns and not columns:
        raise ValueError("SwissProt functional_columns requires an explicit or actual header")
    for meaning, column in functional_columns.items():
        if meaning not in {"go", "ko", "pathway", "cog"} or column not in (columns or []):
            raise ValueError(f"Invalid SwissProt functional column mapping: {meaning}={column}")
    audit.source_info[(species, "SwissProt")] = {
        "header": columns, "comments": comments, "schema_declared": columns is not None,
        "functional_columns": functional_columns,
        "threshold_policy": "No new e-value/identity/coverage threshold or best-hit selection applied",
    }
    return records


def pathway_tokens(value, interpro=False):
    if not is_present(value):
        return set()
    if interpro:
        # The IPR pathway column mixes Reactome, MetaCyc, UniPathway, etc.
        # Only explicit KEGG labels or ko/map-prefixed IDs can become KEGG IDs.
        explicit = re.findall(r"\bKEGG\s*:\s*((?:ko|map)?\d{5})(?!\d)", value, re.I)
        prefixed = re.findall(r"(?<![A-Za-z0-9_])(?:ko|map)\d{5}(?!\d)", value)
        return set(explicit + prefixed)
    return {item for item in re.split(r"[\s,;|]+", value.strip()) if is_present(item)}


def split_cog(value):
    return set(value) & set(COG_NAMES) if is_present(value) else set()


def source_functional(rows, source, swiss_spec=None):
    if source == "InterPro":
        return {key: set(rows.get(key, ())) for key in ("go", "ko", "pathway_tokens", "cog")}
    result = {key: set() for key in ("go", "ko", "pathway_tokens", "cog")}
    for row in rows:
        if source == "eggNOG":
            fields = {"go": row.get("GOs", ""), "ko": row.get("KEGG_ko", ""),
                      "pathway": row.get("KEGG_Pathway", ""), "cog": row.get("COG_category", "")}
        else:
            fields = {meaning: row.get(column, "") for meaning, column in (swiss_spec or {}).get("functional_columns", {}).items()}
        result["go"].update(GO_RE.findall(fields.get("go", "")))
        result["ko"].update(KO_RE.findall(fields.get("ko", "")))
        result["pathway_tokens"].update(pathway_tokens(fields.get("pathway", ""), source == "InterPro"))
        result["cog"].update(split_cog(fields.get("cog", "")))
    return result


def build_gene_records(species, proteins, source_rows, available, go, kegg, swiss_spec, go_audit, pathway_audit, interpro_path=""):
    """Return one record per FASTA protein/gene; source rows remain available raw."""
    result = []
    for query_id, protein in proteins.items():
        egg_rows = source_rows["eggNOG"].get(query_id, [])
        ipr_rows = source_rows["InterPro"].get(query_id, {})
        swiss_rows = source_rows["SwissProt"].get(query_id, [])
        by_source = {
            "eggNOG": source_functional(egg_rows, "eggNOG"),
            "InterPro": source_functional(ipr_rows, "InterPro"),
            "SwissProt": source_functional(swiss_rows, "SwissProt", swiss_spec),
        }
        ipr_ids = ipr_rows.get("ipr_ids", set())
        signatures = ipr_rows.get("signatures", set())
        pfam_ids = ipr_rows.get("pfam_ids", set())
        egg_pfams = {value for row in egg_rows for value in PFAM_RE.findall(row.get("PFAMs", ""))}
        hits = {"eggNOG": any(is_present(row.get("seed_ortholog", "")) for row in egg_rows),
                "InterPro": bool(ipr_ids), "SwissProt": bool(swiss_rows)}
        raw_go = set().union(*(values["go"] for values in by_source.values()))
        canonical_go = set()
        slim_ids = set()
        unresolved_go = set()
        for identifier in sorted(raw_go):
            mapped = go.mapping(identifier)
            if mapped["canonical_id"]:
                canonical_go.add(mapped["canonical_id"])
                slim_ids.update(mapped["slim_ids"])
            else:
                unresolved_go.add(identifier)
            key = (species, identifier)
            if key not in go_audit:
                go_audit[key] = dict(species=species, input_GO_ID=identifier,
                    canonical_GO_ID=mapped["canonical_id"] or "", status=mapped["status"],
                    resolution_chain=";".join(mapped["chain"]), consider=joined(mapped["consider"]),
                    canonical_name=mapped["name"], canonical_namespace=mapped["namespace"],
                    slim_ids=joined(mapped["slim_ids"]), genes=set(), sources=set())
            go_audit[key]["genes"].add(protein["gene_id"])
            go_audit[key]["sources"].update(source for source, values in by_source.items() if identifier in values["go"])
        kos = set().union(*(values["ko"] for values in by_source.values()))
        paths = set()
        paths_by_source = {}
        for source, values in by_source.items():
            normalized_source_paths = set()
            for token in sorted(values["pathway_tokens"]):
                normalized = normalize_pathway(token)
                status = "known" if normalized in kegg.pathways else "unknown_pathway" if normalized else "invalid_pathway_id"
                key = (species, token)
                if key not in pathway_audit:
                    pathway_audit[key] = dict(species=species, input_pathway_ID=token,
                        normalized_pathway_ID=normalized or "", status=status, genes=set(), sources=set())
                pathway_audit[key]["genes"].add(protein["gene_id"])
                pathway_audit[key]["sources"].add(source)
                if normalized:
                    paths.add(normalized)
                    normalized_source_paths.add(normalized)
            paths_by_source[source] = normalized_source_paths
        hierarchy = kegg.hierarchy_for(paths)
        l2 = {(item["A_id"], item["B_id"]) for item in hierarchy}
        cog = set().union(*(values["cog"] for values in by_source.values()))
        row = dict(species=species, **protein,
            eggNOG=flag(hits["eggNOG"], available["eggNOG"]),
            InterProScan_any_signature=flag(signatures, available["InterPro"]),
            InterPro=flag(hits["InterPro"], available["InterPro"]),
            Pfam=flag(pfam_ids, available["InterPro"]),
            SwissProt=flag(hits["SwissProt"], available["SwissProt"]),
            any_upset_source=flag(any(hits[source] for source in SOURCE_ORDER if available[source]), any(available.values())),
            eggNOG_row_count=len(egg_rows) if available["eggNOG"] else "NA",
            InterProScan_row_count=ipr_rows.get("row_count", 0) if available["InterPro"] else "NA",
            SwissProt_row_count=len(swiss_rows) if available["SwissProt"] else "NA",
            GO_raw_union=joined(raw_go), GO_canonical_union=joined(canonical_go), GO_unresolved=joined(unresolved_go),
            GO_slim_union=joined(slim_ids), KEGG_KO_union=joined(kos), KEGG_pathway_union=joined(paths),
            KEGG_pathways_from_KO_reference=joined(kegg.pathways_for_kos(kos)),
            KEGG_L2_IDs=joined(f"{a}:{b}" for a, b in l2), COG_categories=joined(cog),
            InterPro_accessions=joined(ipr_ids), InterProScan_signatures=joined(signatures),
            Pfam_InterProScan_accessions=joined(pfam_ids), Pfam_eggNOG_accessions=joined(egg_pfams),
            InterProScan_analyses=joined(ipr_rows.get("analyses", ())),
            eggNOG_raw_fields_json=json_text(egg_rows) if available["eggNOG"] else "NA",
            InterProScan_domain_records_json=json_text(ipr_rows.get("domains", [])) if available["InterPro"] else "NA",
            InterProScan_raw_row_numbers=json_text(ipr_rows.get("line_numbers", [])) if available["InterPro"] else "NA",
            InterProScan_raw_source_file=interpro_path if available["InterPro"] else "NA",
            SwissProt_raw_fields_json=json_text(swiss_rows) if available["SwissProt"] else "NA",
        )
        for source in SOURCE_ORDER:
            row[f"GO_{source}"] = joined(by_source[source]["go"]) if available[source] else "NA"
            row[f"KO_{source}"] = joined(by_source[source]["ko"]) if available[source] else "NA"
            row[f"KEGG_pathway_{source}"] = joined(paths_by_source[source]) if available[source] else "NA"
        for field in EGG_EXPORT:
            row[f"eggNOG_{field}"] = json_text(sorted({item[field] for item in egg_rows if field in item and is_present(item[field])})) if available["eggNOG"] else "NA"
        row["_sets"] = dict(raw_go=raw_go, canonical_go=canonical_go, slim=slim_ids, ko=kos,
                            pathways=paths, l2=l2, cog=cog, ipr=ipr_ids, signatures=signatures, pfam=pfam_ids)
        result.append(row)
    return result


def count_record(species, source, count, total, available=True):
    return {"species": species, "source": source, "gene_count": count if available else "NA",
            "total_genes": total, "percent_total_genes": percent(count, total) if available else "NA"}


def aggregate_species(species, records, available, go, kegg, source_info, swiss_spec):
    total = len(records)
    order = [source for source in SOURCE_ORDER if available[source]]
    intersections = Counter("".join("1" if row[source] == "true" else "0" for source in order) for row in records)
    upset = []
    for bits in itertools.product("01", repeat=len(order)):
        pattern = "".join(bits)
        count = intersections[pattern]
        upset.append(dict(species=species, pattern=pattern,
            sources=";".join(source for source, bit in zip(order, bits) if bit == "1") or "None",
            gene_count=count, total_genes=total, percent_total_genes=percent(count, total), source_order=";".join(order)))
    assert sum(item["gene_count"] for item in upset) == total
    set_sizes = [count_record(species, source, sum(row[source] == "true" for row in records), total) for source in order]
    coverage = [count_record(species, source, sum(row[source] == "true" for row in records), total, available[source]) for source in SOURCE_ORDER]
    coverage += [count_record(species, source, sum(row[source] == "true" for row in records), total, available["InterPro"])
                 for source in ("InterProScan_any_signature", "Pfam")]
    widths = source_info.get((species, "InterPro"), {}).get("row_width_counts", {})
    ipr_go_available = any(int(width) >= 14 for width in widths)
    ipr_pathway_available = any(int(width) >= 15 for width in widths)
    sw_fields = swiss_spec.get("functional_columns", {})
    go_available = available["eggNOG"] or ipr_go_available or (available["SwissProt"] and "go" in sw_fields)
    ko_available = available["eggNOG"] or (available["SwissProt"] and "ko" in sw_fields)
    path_available = available["eggNOG"] or ipr_pathway_available or (available["SwissProt"] and "pathway" in sw_fields)
    cog_available = available["eggNOG"] or (available["SwissProt"] and "cog" in sw_fields)
    metrics = [("Any_UpSet_source", None, any(available.values())), ("GO_raw", "raw_go", go_available),
               ("GO_resolved", "canonical_go", go_available), ("GO_slim", "slim", go_available),
               ("KEGG_KO", "ko", ko_available), ("KEGG_pathway", "pathways", path_available),
               ("KEGG_L2", "l2", path_available), ("COG", "cog", cog_available)]
    for metric, key, present in metrics:
        count = sum(bool(row["_sets"][key]) if key else row["any_upset_source"] == "true" for row in records)
        coverage.append(count_record(species, metric, count, total, present))
    counts = {key: Counter(value for row in records for value in row["_sets"][key]) for key in ("slim", "pathways", "l2", "cog")}
    slim = [dict(species=species, ontology=ONTOLOGY_LABELS[go.terms[key].namespace], GO_ID=key, GO_term=go.terms[key].name,
                 gene_count=counts["slim"][key] if go_available else "NA", total_genes=total,
                 percent_total_genes=percent(counts["slim"][key], total) if go_available else "NA")
            for key in sorted(go.slim_ids, key=lambda key: (ONTOLOGY_LABELS[go.terms[key].namespace], key))]
    l2_names = {(row["A_id"], row["B_id"]): row for row in kegg.records}
    l2 = [dict(species=species, A_id=a, A_name=l2_names[(a, b)]["A_name"], B_id=b, B_name=l2_names[(a, b)]["B_name"],
               gene_count=counts["l2"][(a, b)] if path_available else "NA", total_genes=total,
               percent_total_genes=percent(counts["l2"][(a, b)], total) if path_available else "NA")
          for a, b in sorted(l2_names)]
    pathways = []
    for pathway_id in sorted(set(kegg.pathways) | set(counts["pathways"])):
        hierarchy = kegg.pathways.get(pathway_id, ())
        pathways.append(dict(species=species, KEGG_pathway_ID=pathway_id,
            KEGG_pathway_name=hierarchy[0]["pathway_name"] if hierarchy else "Unresolved in local hierarchy",
            gene_count=counts["pathways"][pathway_id] if path_available else "NA", total_genes=total,
            percent_total_genes=percent(counts["pathways"][pathway_id], total) if path_available else "NA",
            A_ids=joined(row["A_id"] for row in hierarchy), B_ids=joined(row["B_id"] for row in hierarchy),
            hierarchy_status="known" if hierarchy else "unknown"))
    cog = [dict(species=species, COG_category=key, description=COG_NAMES[key],
                gene_count=counts["cog"][key] if cog_available else "NA", total_genes=total,
                percent_total_genes=percent(counts["cog"][key], total) if cog_available else "NA") for key in sorted(COG_NAMES)]
    lengths = [dict(species=species, gene_id=row["gene_id"], protein_id=row["protein_id"], protein_aa=row["protein_aa"],
                    terminal_stop_count=row["terminal_stop_count"], internal_stop_count=row["internal_stop_count"]) for row in records]
    values = sorted(row["protein_aa"] for row in records)
    length_summary = dict(species=species, total_genes=total, min_aa=min(values), median_aa=statistics.median(values),
                          mean_aa=round(statistics.mean(values), 6), max_aa=max(values), total_aa=sum(values),
                          proteins_with_terminal_stop=sum(row["terminal_stop_count"] > 0 for row in records),
                          proteins_with_internal_stop=sum(row["internal_stop_count"] > 0 for row in records))
    return {"functional_coverage.tsv": coverage, "upset_intersections.tsv": upset, "upset_set_sizes.tsv": set_sizes,
            "goslim_counts.tsv": slim, "kegg_level2_counts.tsv": l2, "kegg_pathway_counts.tsv": pathways,
            "cog_counts.tsv": cog, "protein_lengths.tsv": lengths, "protein_length_summary.tsv": [length_summary]}


def normalize_spec(value, base):
    if value is None or value == "":
        return None
    spec = {"path": value} if isinstance(value, str) else dict(value)
    if not spec.get("path"):
        raise ValueError("Source object requires path")
    path = Path(spec["path"]).expanduser()
    spec["path"] = str((base / path).resolve() if not path.is_absolute() else path.resolve())
    return spec


def export_audits(out, audit, species_data):
    unknown_rows = []
    id_rows = []
    for species, data in species_data.items():
        universe = set(data["proteins"])
        for source in SOURCE_ORDER:
            counts = audit.query_counts[(species, source)]
            unknown = set(counts) - universe
            for query_id in sorted(unknown):
                unknown_rows.append(dict(species=species, source=source, query_id=query_id,
                    input_row_count=counts[query_id], first_line=audit.query_first_line[(species, source, query_id)]))
            present = data["available"][source]
            id_rows.append(dict(species=species, source=source, provided=flag(present), input_rows=sum(counts.values()) if present else "NA",
                unique_query_ids=len(counts) if present else "NA", matched_query_ids=len(set(counts) & universe) if present else "NA",
                unknown_query_ids=len(unknown) if present else "NA", repeated_query_rows=sum(counts.values()) - len(counts) if present else "NA",
                fasta_queries_without_source_rows=len(universe - set(counts)) if present else "NA", total_fasta_genes=len(universe)))
            if unknown:
                audit.issue(species, source, "query_ids_outside_fasta", f"{len(unknown)} distinct query IDs not present in representative FASTA; see unknown_query_ids.tsv")
    write_tsv(out / "query_id_audit.tsv", ["species", "source", "provided", "input_rows", "unique_query_ids", "matched_query_ids", "unknown_query_ids", "repeated_query_rows", "fasta_queries_without_source_rows", "total_fasta_genes"], id_rows)
    write_tsv(out / "unknown_query_ids.tsv", ["species", "source", "query_id", "input_row_count", "first_line"], unknown_rows)
    write_tsv(out / "validation_issues.tsv", ["species", "source", "severity", "kind", "line", "query_id", "detail"], audit.issues)
    write_tsv(out / "interpro_length_audit.tsv", ["species", "comparison", "input_rows"],
              [dict(species=species, comparison=comparison, input_rows=count) for (species, comparison), count in sorted(audit.length_checks.items())])


def run_analysis(manifest_path, out_path):
    manifest_path = Path(manifest_path).resolve()
    config = json.loads(manifest_path.read_text(encoding="utf-8"))
    species_config = config.get("species")
    if not isinstance(species_config, dict) or not species_config:
        raise ValueError("Manifest requires a nonempty species object")
    out = Path(out_path).resolve()
    out.mkdir(parents=True, exist_ok=True)
    existing = [name for name in FUNCTIONAL_RESULTS if (out / name).exists()]
    if existing:
        raise ValueError("Use a new/empty output directory; existing functional results: " + ", ".join(existing))
    base = manifest_path.parent
    go_path = normalize_spec(config.get("go", "go-basic.obo"), base)["path"]
    kegg_path = normalize_spec(config.get("kegg", "ko00001.keg"), base)["path"]
    go, kegg = load_go(go_path), load_kegg(kegg_path)
    audit = Audit()
    species_data = {}
    inputs = [{"role": "manifest", "path": str(manifest_path)}, {"role": "GO", "path": go_path}, {"role": "KEGG", "path": kegg_path}]
    for species, raw_spec in species_config.items():
        if not isinstance(raw_spec, dict) or not raw_spec.get("proteins"):
            raise ValueError(f"Species {species!r} requires proteins")
        specs = {key: normalize_spec(raw_spec.get(key), base) for key in ("proteins", "eggnog", "interpro", "swissprot")}
        source_specs = {source: specs[SOURCE_KEYS[source]] for source in SOURCE_ORDER}
        available = {source: spec is not None for source, spec in source_specs.items()}
        try:
            with open_text(specs["proteins"]["path"]) as handle:
                proteins = read_fasta(handle, species, audit)
        except (OSError, UnicodeError) as exc:
            audit.issue(species, "proteins", "unreadable_file", str(exc))
            proteins = {}
        sources = {source: {} for source in SOURCE_ORDER}
        for source, spec in source_specs.items():
            if spec is None:
                continue
            try:
                with open_text(spec["path"]) as handle:
                    if source == "eggNOG":
                        sources[source] = read_eggnog(handle, species, proteins, audit)
                    elif source == "InterPro":
                        sources[source] = read_interpro(handle, species, proteins, audit)
                    else:
                        sources[source] = read_swissprot(handle, species, proteins, audit, spec)
            except (OSError, UnicodeError, ValueError) as exc:
                audit.issue(species, source, "source_read_error", str(exc))
        species_data[species] = dict(proteins=proteins, sources=sources, available=available,
                                     swiss_spec=specs["swissprot"] or {}, source_specs=source_specs)
        for key, spec in specs.items():
            if spec is not None:
                inputs.append({"species": species, "role": key, "path": spec["path"]})
    export_audits(out, audit, species_data)
    for item in inputs:
        path = Path(item["path"])
        if path.is_file():
            item.update(bytes=path.stat().st_size, sha256=sha256_file(path))
    report = {
        "status": "failed" if audit.errors() else "validated",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": inputs, "script_sha256": sha256_file(__file__),
        "ontology_script_sha256": sha256_file(Path(__file__).with_name("ontology.py")),
        "go_release": go.header.get("data-version", []),
        "source_metadata": {f"{species}/{source}": info for (species, source), info in audit.source_info.items()},
        "source_availability": {species: data["available"] for species, data in species_data.items()},
        "validation_error_count": len(audit.errors()),
        "validation_warning_count": sum(item["severity"] == "warning" for item in audit.issues),
        "methods": {
            "unit": "unique gene= ID; one representative FASTA protein per gene",
            "ID_policy": "exact query identity only; any database query outside FASTA fails",
            "upset_source_order": list(SOURCE_ORDER),
            "eggNOG_hit": "at least one nonmissing seed_ortholog",
            "InterPro_hit": "at least one real IPR accession in standard column 12",
            "InterProScan_any_signature": "nonmissing signature accession in standard column 5",
            "Pfam_hit": "InterProScan analysis=Pfam with nonmissing signature accession",
            "SwissProt_hit": "at least one input row; thresholds and best-hit selection not inferred",
            "GO": "source union, alt/single replaced_by normalization; all is_a/part_of generic-slim ancestors excluding roots; per-gene bins deduplicated",
            "GO_negative_annotations": "Inputs must contain positive annotation fields; GAF NOT-qualified records are not supported",
            "KEGG": "direct annotated pathway union only; ko/map aliases normalized; all pathway A/B memberships retained; each gene counted once per bin",
            "KO_inference": "KO-derived pathways exported separately and excluded from primary pathway/L2 counts",
            "percentage_denominator": "all representative FASTA genes, including unannotated",
            "missing_source": "NA; no-hit in a provided source=false",
            "protein_aa": "sequence length excluding '*' markers; no inferred CDS length is reported",
            "InterPro_length": "exact aa length; +1 tolerated only with one known terminal '*' and no internal stops",
        },
    }
    if audit.errors():
        (out / "analysis_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return 2
    combined = defaultdict(list)
    gene_rows = []
    go_audit, pathway_audit = {}, {}
    for species, data in species_data.items():
        records = build_gene_records(species, data["proteins"], data["sources"], data["available"], go, kegg,
                                     data["swiss_spec"], go_audit, pathway_audit,
                                     (data["source_specs"]["InterPro"] or {}).get("path", ""))
        tables = aggregate_species(species, records, data["available"], go, kegg, audit.source_info, data["swiss_spec"])
        for name, rows in tables.items():
            combined[name].extend(rows)
        for row in records:
            gene_rows.append({key: value for key, value in row.items() if key != "_sets"})
    leading_fields = ["species", "gene_id", "protein_id", "protein_aa", "eggNOG",
                      "InterProScan_any_signature", "InterPro", "Pfam", "SwissProt"]
    gene_fields = leading_fields + [key for key in gene_rows[0] if key not in leading_fields]
    write_tsv(out / "gene_functional_annotations.tsv", gene_fields, gene_rows)
    for name, rows in combined.items():
        if rows:
            write_tsv(out / name, list(rows[0]), rows)
    for filename, items, fields in [
        ("observed_go_audit.tsv", go_audit, ["species", "input_GO_ID", "canonical_GO_ID", "status", "resolution_chain", "consider", "canonical_name", "canonical_namespace", "slim_ids", "gene_count", "sources"]),
        ("observed_pathway_audit.tsv", pathway_audit, ["species", "input_pathway_ID", "normalized_pathway_ID", "status", "gene_count", "sources"]),
    ]:
        rows = []
        for key in sorted(items):
            item = items[key]
            rows.append({**{field: value for field, value in item.items() if field not in {"genes", "sources"}},
                         "gene_count": len(item["genes"]), "sources": joined(item["sources"])})
        write_tsv(out / filename, fields, rows)
    report["status"] = "complete"
    report["representative_genes"] = {species: len(data["proteins"]) for species, data in species_data.items()}
    report["outputs"] = {path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in sorted(out.glob("*.tsv"))}
    (out / "analysis_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


def self_test():
    """In-memory contract checks only. No fake sample files/tables are written."""
    audit = Audit()
    proteins = read_fasta(io.StringIO(">p1 gene=g1\nMAAA*\n>p2 gene=g2\nMAAA\n>p3 gene=g3\nMAAA\n"), "test", audit)
    egg = read_eggnog(io.StringIO(
        "#query\tseed_ortholog\tCOG_category\tGOs\tKEGG_ko\tKEGG_Pathway\n"
        "p1\tx\tK\tGO:0000004\tko:K00001\tko00010,map00010\n"
        "p2\t-\t-\t-\t-\t-\n"), "test", proteins, audit)
    ipr = read_interpro(io.StringIO(
        "p1\tx\t5\tPfam\tPF00001\tdomain\t1\t4\t1e-5\tT\tdate\tIPR000001\tdesc\tGO:0000002|GO:0000002\tKEGG: 00010\n"
        "p1\tx\t5\tPfam\tPF00001\tdomain\t1\t4\t1e-5\tT\tdate\tIPR000001\tdesc\tGO:0000002\tKEGG: 00010\n"
        "p2\tx\t4\tPfam\tPF00002\tdomain\t1\t4\t1e-5\tT\tdate\t-\t-\t-\t-\n"), "test", proteins, audit)
    go = GOOntology({
        "GO:0000001": GOTerm("GO:0000001", "slim", "biological_process", subsets=("goslim_generic",)),
        "GO:0000002": GOTerm("GO:0000002", "child", "biological_process", alt_ids=("GO:0000003",), parents=(("is_a", "GO:0000001"),)),
        "GO:0000004": GOTerm("GO:0000004", "old", "biological_process", is_obsolete=True, replaced_by=("GO:0000003",)),
    }, {}, "memory")
    from ontology import KEGGHierarchy
    kegg = KEGGHierarchy([
        dict(pathway_id="map00010", pathway_name="path", A_id="09100", A_name="Metabolism", B_id="09101", B_name="category1"),
        dict(pathway_id="map00010", pathway_name="path", A_id="09100", A_name="Metabolism", B_id="09102", B_name="category2"),
    ], {"K00001": {"map00010"}}, "memory", 0)
    available = {"eggNOG": True, "InterPro": True, "SwissProt": False}
    rows = build_gene_records("test", proteins, {"eggNOG": egg, "InterPro": ipr, "SwissProt": {}}, available, go, kegg, {}, {}, {})
    tables = aggregate_species("test", rows, available, go, kegg, audit.source_info, {})
    assert not audit.errors(), audit.issues
    assert len(rows) == 3 and rows[0]["GO_canonical_union"] == "GO:0000002"
    assert rows[0]["KEGG_pathway_union"] == "map00010" and rows[0]["SwissProt"] == "NA"
    assert rows[1]["eggNOG"] == "false" and rows[1]["InterProScan_any_signature"] == "true" and rows[1]["InterPro"] == "false"
    assert {row["pattern"]: row["gene_count"] for row in tables["upset_intersections.tsv"]} == {"00": 2, "01": 0, "10": 0, "11": 1}
    assert tables["goslim_counts.tsv"][0]["gene_count"] == 1
    assert [row["gene_count"] for row in tables["kegg_level2_counts.tsv"]] == [1, 1]
    assert audit.length_checks[("test", "terminal_stop_plus_one")] == 2
    invalid_audit = Audit()
    read_fasta(io.StringIO(">p1 gene=g\nMAA\n>p2 gene=g\nMAA\n"), "test", invalid_audit)
    assert any(item["kind"] == "multiple_representatives_per_gene" for item in invalid_audit.errors())
    # Exercise the actual out-of-universe failure branch without writing any
    # synthetic TSV or sample file, even temporarily.
    from unittest.mock import patch
    bounds_audit = Audit()
    bounds_audit.query("test", "eggNOG", "not_in_fasta", 2)
    captured_audits = {}
    def capture_table(path, fields, records):
        captured_audits[path.name] = list(records)
    with patch(__name__ + ".write_tsv", side_effect=capture_table):
        export_audits(Path("."), bounds_audit, {"test": {"proteins": proteins, "available": available}})
    assert any(item["kind"] == "query_ids_outside_fasta" for item in bounds_audit.errors())
    assert captured_audits["unknown_query_ids.tsv"][0]["query_id"] == "not_in_fasta"
    print("In-memory functional-analysis self-test passed; no sample output was written.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", help="JSON input manifest")
    parser.add_argument("--out", help="New output directory")
    parser.add_argument("--self-test", action="store_true", help="Run in-memory contract checks without creating files")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.manifest or not args.out:
        parser.error("--manifest and --out are required unless --self-test is used")
    try:
        status = run_analysis(args.manifest, args.out)
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if status:
        print("Validation failed; inspect validation_issues.tsv and unknown_query_ids.tsv. No functional tables were generated.", file=sys.stderr)
    else:
        print(f"Functional analysis complete: {Path(args.out).resolve()}")
    return status


if __name__ == "__main__":
    sys.exit(main())
