import pandas as pd
import matplotlib.pyplot as plt
import os
from glob import glob
import argparse

## Parser
parser = argparse.ArgumentParser(
    description="Process bracken files and plot genus-level composition"
)

parser.add_argument(
    "-i",
    "--input_dir",
    required=True,
    help="Directory containing *_G.txt files"
)

parser.add_argument(
    "-o",
    "--output_dir",
    required=True,
    help="Directory where plots will be saved"
)

args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

## -------------------------
## Load and preprocess files
## -------------------------

all_dfs = []

for filepath in glob(os.path.join(args.input_dir, "*_G.txt")):

    fname = os.path.basename(filepath)

    sample = fname.replace("_G.txt", "")
    sample = sample.replace("trimmed_", "")

    df = pd.read_csv(filepath, sep="\t")

    # merge <1% into Others
    low_abundance = df["fraction_total_reads"] < 0.01

    others_fraction = df.loc[
        low_abundance,
        "fraction_total_reads"
    ].sum()

    df = df.loc[
        ~low_abundance,
        ["name", "fraction_total_reads"]
    ].copy()

    if others_fraction > 0:
        df.loc[len(df)] = ["Others", others_fraction]

    df["Sample"] = sample

    all_dfs.append(df)

## -------------------------
## Build matrix
## -------------------------

combined_df = pd.concat(all_dfs, ignore_index=True)

plot_df = combined_df.pivot_table(
    index="Sample",
    columns="name",
    values="fraction_total_reads",
    aggfunc="sum",
    fill_value=0
)

# normalize per sample
plot_df = plot_df.div(plot_df.sum(axis=1), axis=0)

## -------------------------
## Collapse rare genera if >11
## -------------------------

n_genera = len([c for c in plot_df.columns if c != "Others"])

if n_genera > 11:

    mean_abundance = plot_df.mean(axis=0)

    if "Others" in mean_abundance.index:
        mean_abundance = mean_abundance.drop("Others")

    top_genera = (
        mean_abundance
        .sort_values(ascending=False)
        .head(11)
        .index
        .tolist()
    )

    low_genera = [
        g for g in plot_df.columns
        if g not in top_genera and g != "Others"
    ]

    if "Others" not in plot_df.columns:
        plot_df["Others"] = 0

    plot_df["Others"] += plot_df[low_genera].sum(axis=1)

    plot_df = plot_df[top_genera + ["Others"]]

## -------------------------
## Final ordering
## -------------------------

plot_df = plot_df.sort_index()

cols = [c for c in plot_df.columns if c != "Others"] + ["Others"]
plot_df = plot_df[cols]

## -------------------------
## Colors (Set3 only)
## -------------------------
Set3 = [
    "#8DD3C7",
    "#FFFFB3",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#FCCDE5",
    "#BC80BD",
    "#CCEBC5",
    "#FFED6F",
    "#D9D9D9"
]


colors = []

i = 0
for genus in plot_df.columns:

    if genus == "Others":
        colors.append("#D9D9D9")
    else:
        colors.append(Set3[i])
        i += 1

## -------------------------
## Plot
## -------------------------

ax = plot_df.plot(
    kind="barh",
    stacked=True,
    figsize=(14, 10),
    color=colors
)

fig = ax.get_figure()

plt.xlabel("Relative Abundance")
plt.ylabel("Sample")
plt.title("Normalized Microbial Community Composition (Genus Level)")

plt.tight_layout(rect=[0, 0, 0.8, 1])

plt.legend(
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    title="Genus",
    frameon=False
)

## -------------------------
## Save
## -------------------------

for ext in ["png"]:
    fig.savefig(
        os.path.join(args.output_dir, f"stacked_barplot.{ext}"),
        dpi=300,
        bbox_inches="tight"
    )

print(f"Plots saved to {args.output_dir}")