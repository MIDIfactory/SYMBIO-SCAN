#!/usr/bin/env python3

from collections import defaultdict
import argparse
import os
import sys

# ----------------------------
# Argument parsing
parser = argparse.ArgumentParser(description="Rename final OTU consensus headers using groups.tsv and final taxonomy")
parser.add_argument("-i", "--input_file", required=True, help="Final OTUs multifasta")
parser.add_argument("-t", "--final_tax", required=True, help="Taxonomy assigned to the final clusters")
parser.add_argument("-g", "--groups_file", required=True, help="groups.tsv file")
parser.add_argument("-s", "--sample", required=True, help="Sample name")
parser.add_argument("-o", "--outdir", required=True, help="Output directory")
parser.add_argument("-c", "--collapse_genus", default="False", help="Collapse clusters with same genus (True/False)"
)
args = parser.parse_args()

fasta_file = args.input_file
final_tax_file = args.final_tax
groups_file = args.groups_file
sample = args.sample
outdir = args.outdir
collapse_genus = args.collapse_genus.lower() == "true"

os.makedirs(outdir, exist_ok=True)

# ----------------------------
# 1) Parse groups.tsv (reference for cluster codes)
clusters = defaultdict(lambda: {"otus": [], "size": 0, "final_taxonomy": None, "ref_identity": None})

with open(groups_file) as f:
    header = next(f)  # skip header
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue
        cluster_code = fields[2]  # numeric cluster code
        otu = fields[5]
        size = int(fields[6])

        clusters[cluster_code]["otus"].append(otu)
        clusters[cluster_code]["size"] += size  # cumulative size

# ----------------------------
# 2) Parse final taxonomy file
# Map taxonomy using numeric cluster code extracted from first field
with open(final_tax_file) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        cluster_full = fields[0]                # e.g., Ouky169_cluster_0
        cluster_code = cluster_full.split("_")[-1]  # extract numeric code
        taxonomy = fields[2]
        ref_identity = fields[3]
        if cluster_code in clusters:
            clusters[cluster_code]["final_taxonomy"] = taxonomy
            clusters[cluster_code]["ref_identity"] = ref_identity
        else:
            print(f"WARNING: cluster {cluster_full} (code {cluster_code}) in final taxonomy not found in groups.tsv", file=sys.stderr)


# ----------------------------
# 3) Resolve duplicate taxonomies
taxonomy_to_cluster = {}
discarded_clusters = set()

for cluster_code, data in clusters.items():
    tax = data["final_taxonomy"]
    if not tax:
        continue  # skip clusters without taxonomy
    if tax not in taxonomy_to_cluster:
        taxonomy_to_cluster[tax] = cluster_code
    else:
        existing = taxonomy_to_cluster[tax]
        if clusters[cluster_code]["size"] > clusters[existing]["size"]:
            discarded_clusters.add(existing)
            taxonomy_to_cluster[tax] = cluster_code
        else:
            discarded_clusters.add(cluster_code)


# ----------------------------
# 3b) Optional: resolve duplicate genera

def extract_genus(taxonomy):
    """
    Extract genus based on:
    - find family ending with 'eae'
    - genus is the next taxonomic rank
    """
    if not taxonomy:
        return None

    parts = [p.strip() for p in taxonomy.split(";")]

    for i, level in enumerate(parts):
        if level.endswith("eae") and i + 1 < len(parts):
            return parts[i + 1]  # genus
    return None


if collapse_genus:
    genus_to_cluster = {}

    for cluster_code, data in clusters.items():
        if cluster_code in discarded_clusters:
            continue  # already discarded

        genus = extract_genus(data["final_taxonomy"])
        if not genus:
            continue

        if genus not in genus_to_cluster:
            genus_to_cluster[genus] = cluster_code
        else:
            existing = genus_to_cluster[genus]
            if clusters[cluster_code]["size"] > clusters[existing]["size"]:
                discarded_clusters.add(existing)
                genus_to_cluster[genus] = cluster_code
            else:
                discarded_clusters.add(cluster_code)

# ----------------------------
# 4) Read final OTUs FASTA
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

# ----------------------------
# 5) Write renamed multifasta
out_fasta = os.path.join(outdir, f"{sample}_high_quality_seqs.fasta")
discarded_out_fasta = os.path.join(outdir, f"{sample}_sus_sequences.fasta")

with open(out_fasta, "w") as out_kept, open(discarded_out_fasta, "w") as out_discarded:
    for header, sequence in fasta.items():
        # Extract numeric cluster code from FASTA header
        if "cluster_" not in header:
            print(f"WARNING: cannot parse cluster from header {header}", file=sys.stderr)
            continue
        cluster_code = header.split("cluster_")[-1]

        if cluster_code not in clusters:
            print(f"WARNING: cluster {cluster_code} not found in groups.tsv", file=sys.stderr)
            continue

        otus = ";".join(clusters[cluster_code]["otus"])
        tot_size = clusters[cluster_code]["size"]
        final_taxonomy = clusters[cluster_code]["final_taxonomy"]
        ref_identity = clusters[cluster_code]["ref_identity"]
        # Build consistent cluster name using sample prefix
        fasta_cluster_name = f"{sample}_cluster_{cluster_code}"

        new_header = (
            f">{sample}\t{fasta_cluster_name}\t"
            f"{otus}\t"
            f"size={tot_size}\t"
            f"{ref_identity}\t"
            f"putative {final_taxonomy}"
        )

        if cluster_code in discarded_clusters:
            out_discarded.write(new_header + "\n")
            out_discarded.write(sequence + "\n")
        else:
            out_kept.write(new_header + "\n")
            out_kept.write(sequence + "\n")

print(f"Renamed final OTUs written to {out_fasta}")
print(f"Sus sequences written to {discarded_out_fasta}")
