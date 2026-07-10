#!/usr/bin/env python3

import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(
    description="Create an OTU-like table from summary fasta TSV file"
)

parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Input summary TSV file"
)

parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="Output OTU table TSV file"
)

args = parser.parse_args()

INPUT = args.input
OUTPUT = args.output

# (sample, column_name) -> (size, sequence)
records = {}

samples = set()

# Will store tuples:
# (sort_key, column_name)
taxonomy_info = []

with open(INPUT, "r") as f:
    while True:
        header = f.readline()
        if not header:
            break

        seq = f.readline().strip()
        header = header.strip()

        if not header.startswith(">"):
            raise ValueError(f"Expected header line beginning with '>', got:\n{header}")

        fields = header[1:].split("\t")

        if len(fields) < 6:
            raise ValueError(f"Header has fewer than 6 fields:\n{header}")

        sample = fields[0]

        # size
        m = re.search(r"size=(\d+)", fields[3])
        if m is None:
            raise ValueError(f"Cannot parse size from {fields[3]}")
        size = int(m.group(1))

        # identity
        identity = fields[4]

        # taxonomy field
        full_taxonomy = fields[5]

        if not full_taxonomy.startswith("putative "):
            raise ValueError(f"Unexpected taxonomy field:\n{full_taxonomy}")

        full_taxonomy = full_taxonomy[len("putative "):]

        parts = full_taxonomy.split(" ", 1)

        accession = parts[0]

        # -------------------------------------------------------
        # CASE 1: accession only
        # -------------------------------------------------------
        if len(parts) == 1:

            column = f"{accession} ({identity}%)"

            # sort after all taxa
            sort_key = (1, accession)

        # -------------------------------------------------------
        # CASE 2: accession + taxonomy
        # -------------------------------------------------------
        else:

            taxonomy = parts[1]

            taxa = taxonomy.split(";")

            if len(taxa) < 2:
                raise ValueError(f"Unexpected taxonomy:\n{taxonomy}")

            genus = taxa[-2]
            species = taxa[-1]

            column = f"{genus};{species} {accession} ({identity}%)"

            # sort by genus, species, accession
            sort_key = (0, genus.lower(), species.lower(), accession)

        samples.add(sample)

        taxonomy_info.append((sort_key, column))

        key = (sample, column)

        # keep largest size if duplicated
        if key not in records or size > records[key][0]:
            records[key] = (size, seq)

# -------------------------------------------------------
# Build ordered list of columns
# -------------------------------------------------------

taxonomy_info = sorted(set(taxonomy_info))

columns = [x[1] for x in taxonomy_info]

samples = sorted(samples)

# -------------------------------------------------------
# Build dataframe
# -------------------------------------------------------

df = pd.DataFrame(index=samples, columns=columns)

for (sample, column), (_, seq) in records.items():
    df.loc[sample, column] = seq

df = df.fillna("NA")

df.index.name = "Sample"

df.to_csv(OUTPUT, sep="\t")

print(f"Wrote {OUTPUT}")
print(f"{len(samples)} samples")
print(f"{len(columns)} taxa")
