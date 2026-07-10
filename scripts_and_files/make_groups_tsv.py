#!/usr/bin/env python3
import argparse
import os
import sys


# Argument parsing
parser = argparse.ArgumentParser(description="Create groups.tsv from vsearch UC and BLAST taxonomy")
parser.add_argument("-u", "--uc_file", required=True, help="vsearch .uc file")
parser.add_argument("-t", "--tax_file", required=True, help="BLAST taxonomy file")
parser.add_argument("-s", "--sample", required=True, help="Sample name")
parser.add_argument("-o", "--outdir", required=True, help="Output directory")
args = parser.parse_args()

uc_file = args.uc_file
tax_file = args.tax_file
sample = args.sample
outdir = args.outdir

os.makedirs(outdir, exist_ok=True)
out_tsv = os.path.join(outdir, f"{sample}_groups.tsv")

# 1) Read taxonomy file → OTU → taxonomy
taxonomy = {}

with open(tax_file) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        # OTU id without ;size=X
        otu = fields[0].split(";size=")[0]

        # Tax info - taxonomy and iIdentity with the reference
        taxonomy[otu] = {
        "taxonomy": fields[2],
        "ref_perc_id": fields[3]
        }


# 2) Parse UC file and write groups.tsv, sorted by cluster
rows = []

with open(uc_file) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        record_type = fields[0]
        if record_type not in ("S", "H"):
            continue

        cluster = fields[1]
        length = fields[2]
        perc_id = fields[3]

        # UC label (contains OTU and size)
        label = fields[8]
        if ";size=" not in label:
            continue
        otu, size = label.split(";size=")

        info = taxonomy.get(otu, {"taxonomy": "NA", "ref_perc_id": "NA"})
        tax = info["taxonomy"]
        ref_perc_id = info["ref_perc_id"]

        rows.append([record_type, cluster, length, perc_id, otu, size, ref_perc_id, tax])

# Sort rows by cluster (numerically)
rows.sort(key=lambda x: int(x[1]))

# Write output
with open(out_tsv, "w") as out:
    out.write("sample\tS/H\tcluster\tlength\tperc_id\tOTU\tsize\tref_perc_id\ttaxonomy\n")
    for row in rows:
        out.write("\t".join([sample] + row) + "\n")

print(f"groups.tsv written to {out_tsv} ({len(rows)} rows)")
