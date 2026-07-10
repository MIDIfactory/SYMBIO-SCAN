#!/usr/bin/env python3

from collections import defaultdict
import argparse
import glob
import os
import sys

# Argument parsing
parser = argparse.ArgumentParser(description="Split OTUs into per-cluster FASTA files")
parser.add_argument("-i", "--input_file", required=True, help="Input name")
parser.add_argument("-g", "--groups_file", required=True, help="Groups file")
parser.add_argument("-s", "--sample", required=True, help="Sample name")
parser.add_argument("-o", "--outdir", required=True, help="Output directory for cluster FASTAs")
args = parser.parse_args()

fasta_file = args.input_file
table_file = args.groups_file
sample = args.sample
outdir = args.outdir

os.makedirs(outdir, exist_ok=True)

# 1) Read groups.tsv → cluster → OTU;size
clusters = defaultdict(list)

with open(table_file) as f:
    header = next(f)
    for line in f:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 8:
            continue

        cluster = fields[2]
        otu = fields[5]
        size = fields[6]

        clusters[cluster].append(f"{otu};size={size}")


# 2) Read FASTA into memory
fasta = {}
with open(fasta_file) as f:
    header = None
    seq = []
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                fasta[header] = "".join(seq)
            header = line[1:]
            seq = []
        else:
            seq.append(line)
    if header:
        fasta[header] = "".join(seq)


# 3) Write one FASTA per cluster
for cluster, otus in clusters.items():
    outname = os.path.join(outdir, f"{sample}_cluster_{cluster}.fasta")
    with open(outname, "w") as out:
        for otu in otus:
            if otu not in fasta:
                print(f"WARNING: {otu} not found in FASTA", file=sys.stderr)
                continue
            out.write(f">{otu}\n")
            out.write(fasta[otu] + "\n")

    print(f"Wrote {outname}")
