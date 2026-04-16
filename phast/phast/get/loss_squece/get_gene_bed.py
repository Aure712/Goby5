#!/usr/bin/env python3
from collections import defaultdict

# --- user config ---
input_file = "/fast3/group_crf/home/g21shaoy23/goby5/annotations/cds_phase_fix2/Tra.prechange_longest_fixed.gff3"
output_file = "gene_noUTR.bed"
use_bed_coords = True  # True: BED 0-based start; False: keep GFF 1-based start
# -------------------

def parse_attrs(attr_str):
    attrs = {}
    for item in attr_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs

records = []
feature_type_by_id = {}
parent_by_id = {}
name_by_id = {}

with open(input_file, "r", encoding="utf-8") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attr_str = parts
        attrs = parse_attrs(attr_str)
        rec_id = attrs.get("ID")
        parents = attrs.get("Parent", "")
        parent_ids = [p for p in parents.split(",") if p]

        if rec_id:
            feature_type_by_id[rec_id] = ftype
            parent_by_id[rec_id] = parent_ids
            name_by_id[rec_id] = (
                attrs.get("Name")
                or attrs.get("gene")
                or attrs.get("gene_name")
                or rec_id
            )

        records.append({
            "seqid": seqid,
            "type": ftype,
            "start": int(start),
            "end": int(end),
            "attrs": attrs,
            "parent_ids": parent_ids,
        })

def find_gene_id(feature_id):
    seen = set()
    stack = [feature_id]
    while stack:
        cur = stack.pop()
        if cur in seen:
            continue
        seen.add(cur)
        ftype = feature_type_by_id.get(cur)
        if ftype == "gene":
            return cur
        for p in parent_by_id.get(cur, []):
            stack.append(p)
    return None

gene_ranges = {}
for rec in records:
    if rec["type"] != "CDS":
        continue

    gene_id = None
    for pid in rec["parent_ids"]:
        gene_id = find_gene_id(pid)
        if gene_id:
            break

    if gene_id:
        gene_name = name_by_id.get(gene_id, gene_id)
        gene_key = gene_id
    else:
        gene_name = (
            rec["attrs"].get("gene")
            or rec["attrs"].get("gene_name")
            or rec["attrs"].get("Name")
            or (rec["parent_ids"][0] if rec["parent_ids"] else "unknown")
        )
        gene_key = gene_name

    entry = gene_ranges.get(gene_key)
    if not entry:
        gene_ranges[gene_key] = {
            "chrom": rec["seqid"],
            "start": rec["start"],
            "end": rec["end"],
            "name": gene_name,
        }
    else:
        entry["start"] = min(entry["start"], rec["start"])
        entry["end"] = max(entry["end"], rec["end"])

# write BED
with open(output_file, "w", encoding="utf-8") as out:
    for gene_key, info in sorted(gene_ranges.items(), key=lambda x: (x[1]["chrom"], x[1]["start"])):
        start = info["start"] - 1 if use_bed_coords else info["start"]
        out.write(f"{info['chrom']}\t{start}\t{info['end']}\t{info['name']}\n")

print(f"Wrote: {output_file}")
