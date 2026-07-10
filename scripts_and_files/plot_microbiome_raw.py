import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
from glob import glob
import argparse

# -----------------------
# Parse command-line args
# -----------------------
parser = argparse.ArgumentParser(description="Process kraken_raw krona files and plot top 30 genera")
parser.add_argument("-i", "--input_dir", required=True, help="Directory containing krona files")
parser.add_argument("-o", "--output_dir", required=True, help="Directory where plots will be saved")
args = parser.parse_args()

# -----------------------
# Prepare temp folder
# -----------------------
temp_dir = os.path.join(args.output_dir, "temp")
os.makedirs(temp_dir, exist_ok=True)

# -----------------------
# Step 1: Process krona files into genus\tcount
# -----------------------
for filepath in glob(os.path.join(args.input_dir, "*.krona")):
    output_file = os.path.join(temp_dir, os.path.basename(filepath) + ".out")
    with open(filepath, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith("0") or not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) < 2:
                continue
            count = fields[0]
            genus = fields[-1]
            fout.write(f"{genus}\t{count}\n")

# -----------------------
# Step 2: Rename cleaned files
# -----------------------
for fname in os.listdir(temp_dir):
    base = fname.replace(".fastq_bracken_genuses.report.krona.out", "")
    base = base.replace("trimmed_", "")
    old_path = os.path.join(temp_dir, fname)
    new_path = os.path.join(temp_dir, base)
    os.rename(old_path, new_path)

# -----------------------
# Step 3: Load sample files
# -----------------------
file_paths = glob(os.path.join(temp_dir, "*"))
data = {}
for path in file_paths:
    sample_name = os.path.basename(path)
    df = pd.read_csv(path, sep="\t", header=None, names=["genus", sample_name])
    df["genus"] = df["genus"].str.replace("g__", "", regex=False)
    df = df.groupby("genus", as_index=False).sum()
    data[sample_name] = df.set_index("genus")[sample_name]

# -----------------------
# Step 4: Merge and normalize
# -----------------------
merged_df = pd.DataFrame(data).fillna(0).astype(int)
plot_df = merged_df.T
plot_df_normalized = plot_df.div(plot_df.sum(axis=1), axis=0)

# -----------------------
# Step 5: Top 30 genera + "Other"
# -----------------------
top_30_genera = plot_df_normalized.sum().sort_values(ascending=False).head(30).index
others = plot_df_normalized.drop(columns=top_30_genera)
plot_df_top = plot_df_normalized[top_30_genera].copy()
plot_df_top["Other"] = others.sum(axis=1)

# -----------------------
# Step 6: Sort samples alphabetically
# -----------------------
plot_df_top = plot_df_top.sort_index()

# -----------------------
# Step 7: Custom color palette (30 + "Other")
# -----------------------
custom_colors = [
    "#e59d41", "#6069d8", "#5fbd47", "#a15acb", "#bcba39", "#cf51b5",
    "#45c584", "#d74788", "#418533", "#dd4461", "#4fc2b3", "#c4412f",
    "#36bddc", "#dc6e35", "#4880ce", "#839831", "#a594df", "#a98027",
    "#655b9e", "#8bbb6d", "#994d88", "#54a47b", "#dc88be", "#2b7650",
    "#a34758", "#589ed3", "#9b5e30", "#c0a968", "#df8578", "#6b6e2b",
    "#999999"  # Grey for "Other"
]

# -----------------------
# Step 8: Plot
# -----------------------
fig, ax = plt.subplots(figsize=(16, 12))
plot_df_top.plot(kind="barh", stacked=True, ax=ax, color=custom_colors)

ax.set_xlabel("Relative Abundance (0–1)")
ax.set_ylabel("Sample")
ax.set_title("Normalized Microbial Community Composition (Top 30 taxa)")

# Adjust layout and legend
plt.tight_layout(rect=[0, 0, 0.78, 1])
ax.legend(
    bbox_to_anchor=(1.01, 1),
    loc='upper left',
    title="Taxa",
    frameon=False,
    fontsize="small",
    title_fontsize="medium"
)

# -----------------------
# Step 9: Save to output_dir
# -----------------------
plot_path = os.path.join(args.output_dir, "stacked_barplot_top30.pdf")
fig.savefig(plot_path, dpi=300, bbox_inches="tight")
print(f"Plot saved to {plot_path}")

# -----------------------
# Step 10: Cleanup
# -----------------------
shutil.rmtree(temp_dir)
