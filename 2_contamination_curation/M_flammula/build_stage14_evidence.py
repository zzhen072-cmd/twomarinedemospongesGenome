#!/usr/bin/env python3
import csv
import json
import math
import statistics
from pathlib import Path

BASE = Path("/work1/sponge/lan_reassembly_compare")
AD = BASE / "10_fcs_gx/lan_fcs_adaptor.fasta"
GX = BASE / "10_fcs_gx/lan_fcs_gx.fasta"
TIARA = BASE / "13_tiara/lan.tiara.tsv"
KRAKEN = BASE / "11_kraken2/lan.kraken2.out"
BLOB = BASE / "12_blobtoolkit/lan_fcs_adaptor.blobdir"
DB = Path("/public/home/sponge/sci_data/S07/software/standard16")
OUT = BASE / "14_manual_decontam"
OUT.mkdir(parents=True, exist_ok=True)


def norm(value):
    value = value.split()[0]
    return value[4:] if value.startswith("lcl|") else value


def fasta(path):
    records, head, seq = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if head is not None:
                    records[norm(head)] = (head, "".join(seq))
                head, seq = line[1:], []
            else:
                seq.append(line)
        if head is not None:
            records[norm(head)] = (head, "".join(seq))
    return records


records = fasta(AD)
gx_ids = set(fasta(GX))
tiara = {}
with open(TIARA) as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        tiara[norm(row["sequence_id"])] = row

parent, taxname = {}, {}
with open(DB / "nodes.dmp") as handle:
    for line in handle:
        fields = line.split("|")
        parent[int(fields[0])] = int(fields[1])
with open(DB / "names.dmp") as handle:
    for line in handle:
        fields = line.split("|")
        if fields[3].strip() == "scientific name":
            taxname[int(fields[0])] = fields[1].strip()


def domain(taxid):
    seen = set()
    while taxid and taxid not in seen:
        seen.add(taxid)
        if taxid == 2:
            return "bacteria"
        if taxid == 2157:
            return "archaea"
        if taxid == 2759:
            return "eukarya"
        if taxid == 10239:
            return "viruses"
        next_taxid = parent.get(taxid, taxid)
        if next_taxid == taxid:
            break
        taxid = next_taxid
    return "other"


kraken = {}
with open(KRAKEN) as handle:
    for line in handle:
        fields = line.rstrip().split("\t")
        if len(fields) >= 4:
            taxid = int(fields[2]) if fields[2].isdigit() else 0
            kraken[norm(fields[1])] = (
                fields[0],
                taxid,
                taxname.get(taxid, ""),
                domain(taxid),
            )


def vector(name):
    payload = json.load(open(BLOB / name))
    if isinstance(payload, dict):
        for key in ("values", "data"):
            if key in payload:
                return payload[key]
    if isinstance(payload, list):
        return payload
    raise ValueError(name)


ids = [norm(str(value)) for value in vector("identifiers.json")]
lengths = vector("length.json")
gc_values = vector("gc.json")
coverage = vector("lan.ont_to_fcs_adaptor_cov.json")
if not (len(ids) == len(lengths) == len(gc_values) == len(coverage) == len(records)):
    raise ValueError("BlobToolKit vectors and FASTA counts differ")
blob = {
    ids[index]: (int(lengths[index]), float(gc_values[index]), float(coverage[index]))
    for index in range(len(ids))
}

host_cov = [
    value[2]
    for key, value in blob.items()
    if tiara.get(key, {}).get("class_fst_stage") == "eukarya"
    and key in gx_ids
    and value[2] > 0
]
host_gc = [
    value[1]
    for key, value in blob.items()
    if tiara.get(key, {}).get("class_fst_stage") == "eukarya" and key in gx_ids
]
med_cov = statistics.median(host_cov)
med_gc = statistics.median(host_gc)
mad_cov = statistics.median(abs(math.log2(value / med_cov)) for value in host_cov)
mad_gc = statistics.median(abs(value - med_gc) for value in host_gc)

rows, delete, review = [], [], []
for name, (header, seq) in records.items():
    length, gc, cov = blob[name]
    tiara_row = tiara.get(name, {})
    tiara_1 = tiara_row.get("class_fst_stage", "missing")
    tiara_2 = tiara_row.get("class_snd_stage", "n/a")
    k_status, k_taxid, k_name, k_domain = kraken.get(
        name, ("U", 0, "", "unclassified")
    )
    fcs = name not in gx_ids
    pro_tiara = tiara_1 in {"bacteria", "archaea", "prokarya"}
    pro_kraken = k_domain in {"bacteria", "archaea", "viruses"}
    cov_outlier = cov > 0 and abs(math.log2(cov / med_cov)) > max(2.5, 6 * mad_cov)
    gc_outlier = abs(gc - med_gc) > max(0.12 if med_gc <= 1 else 12, 6 * mad_gc)
    support = sum((fcs, pro_tiara, pro_kraken, cov_outlier, gc_outlier))
    conflict = tiara_1 == "eukarya" and not pro_kraken
    if (
        (
            (fcs and (pro_tiara or pro_kraken))
            or (pro_tiara and pro_kraken)
            or (pro_tiara and cov_outlier and gc_outlier)
            or (pro_kraken and cov_outlier and gc_outlier)
        )
        and not (conflict and support < 3)
    ):
        decision = "DELETE"
    elif fcs or pro_tiara or pro_kraken or cov_outlier or gc_outlier:
        decision = "REVIEW"
    else:
        decision = "KEEP"
    reason = ";".join(
        key
        for key, value in (
            ("fcs_gx_removed", fcs),
            ("tiara_prok", pro_tiara),
            ("kraken_prok_or_virus", pro_kraken),
            ("coverage_outlier", cov_outlier),
            ("gc_outlier", gc_outlier),
        )
        if value
    ) or "no_contamination_signal"
    if decision == "DELETE":
        delete.append(name)
    elif decision == "REVIEW":
        review.append(name)
    rows.append(
        [
            name,
            length,
            gc,
            cov,
            fcs,
            tiara_1,
            tiara_2,
            k_status,
            k_taxid,
            k_name,
            k_domain,
            cov_outlier,
            gc_outlier,
            support,
            decision,
            reason,
        ]
    )

header = [
    "contig",
    "length",
    "gc",
    "ont_cov",
    "fcs_gx_removed",
    "tiara_stage1",
    "tiara_stage2",
    "kraken_status",
    "kraken_taxid",
    "kraken_name",
    "kraken_domain",
    "coverage_outlier",
    "gc_outlier",
    "support_count",
    "decision",
    "reason",
]
with open(OUT / "lan_decontam_evidence.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(header)
    writer.writerows(rows)
for filename, values in (
    ("delete_candidates.txt", delete),
    ("review_candidates.txt", review),
):
    (OUT / filename).write_text("".join(value + "\n" for value in values))
with open(OUT / "audit_summary.txt", "w") as handle:
    handle.write(
        f"input_contigs\t{len(records)}\n"
        f"input_bases\t{sum(len(seq) for _, seq in records.values())}\n"
    )
    handle.write(
        f"fcs_gx_removed\t{len(set(records) - gx_ids)}\n"
        f"delete_candidates\t{len(delete)}\n"
        f"review_candidates\t{len(review)}\n"
    )
    handle.write(
        f"delete_bases\t{sum(len(records[name][1]) for name in delete)}\n"
        f"review_bases\t{sum(len(records[name][1]) for name in review)}\n"
    )
    handle.write(f"host_median_cov\t{med_cov}\nhost_median_gc\t{med_gc}\n")

delete_set = set(delete)
with open(OUT / "lan_manual_decontaminated.fasta", "w") as keep, open(
    OUT / "lan_removed_contaminants.fasta", "w"
) as removed:
    for name, (fasta_header, seq) in records.items():
        output = removed if name in delete_set else keep
        output.write(">" + fasta_header + "\n")
        for start in range(0, len(seq), 80):
            output.write(seq[start : start + 80] + "\n")
print((OUT / "audit_summary.txt").read_text(), end="")
