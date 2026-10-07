#!/bin/bash

USAGE="
Usage: ./6.phylogeny.sh -i <data_path> [options]

Required:
  -i  Input directory containing FASTQ files

Optional:
  -o  Output directory (default: basename of input path)
  -t  Number of threads (default: 12)

  -h  Show this help message
"

while getopts "i:o:t:h" OPTION; do
  case $OPTION in
        i) data_path=$OPTARG;;
        o) results=$OPTARG;;
        t) threads=$OPTARG;;
        h) echo "$USAGE"; exit 0;;
        :) printf "Missing argument for option -$OPTARG\n" >&2; exit 1;;
        \?) printf "Invalid option -$OPTARG\n" >&2; exit 1;;
  esac
done

# Validate input
MISSING="is missing but required. Exiting."
if [ -z "${data_path+x}" ]; then
    echo "-i $MISSING"
    echo "$USAGE"
    echo "Enjoy"
    exit 1
fi

# Set defaults
base_data_path=$(basename $data_path)
[ -z "${results+x}" ] && results="${base_data_path}_SymbiontScreener_results"
results_name=$(basename "$results")
[ -z "${threads+x}" ] && threads=12

# Create output directory
mkdir -p "$results"

# Log
mkdir -p "$results/logs"
exec > >(tee -a "$results/logs/6.phylogeny.txt") 2>&1

#:::::::::::::::::::::::::::::::::::::::::::::::::::
# Announce configuration
echo "This is 6.phylogeny.sh"
echo "==== PRELIMIANRY PHYLOGENY ===="
echo "Start time: $(date)"
echo "Input data path: $data_path"
echo "Output directory: $results/6.phylogeny"
echo "Threads: $threads"
echo ""

# Welcoming ducks
cat <<'EOF'
 __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<
 \___)   \___)   \___)   \___)   \___)   \___)   \___)   \___)
quack quack quack quack quack quack quack quack quack quack quack
EOF
echo ""

#:::::::::::::::::::::::::::::::::::::::::::::: 
# Phylogeny out_dir
phylo_out="$results/6.phylogeny"
mkdir -p "$phylo_out"

# Define databases
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SILVA_DB="$SCRIPT_DIR/16S_DB/SILVA138/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta"

# Process the samples all together
for file in "$data_path"/*.fastq; do
  sample=$(basename "$file" .fastq)

  echo "Process sample ${sample}"

  # Merge the input from different samples in a multifasta
  cat \
  "$results/5.amplicons_seqs/$sample/14.final_sequences/${sample}_high_quality_seqs.fasta" | \
  awk -F'\t' '/^>/{print ">"$2; next} {print}' \
  >> "$phylo_out"/"${results_name}_high_quality_sequences.fasta"


  # Extraxt reference sequences
  cut -f2 "$results/5.amplicons_seqs/$sample/12.cons_correction/3.final_otus/${sample}_final_otus_tax.txt" | \
  while IFS= read -r i; do
      samtools faidx "$SILVA_DB" "$i" \
        >> "$phylo_out"/"${results_name}_reference_sequences.fasta"
  done
  
  # Retrieve the names of the reference sequences
  cut -f2,3 "$results/5.amplicons_seqs/$sample/12.cons_correction/3.final_otus/${sample}_final_otus_tax.txt" | sort -u | \
  while read -r A B; do
      sed -i "s/$A/$B/g" "$phylo_out"/"${results_name}_reference_sequences.fasta"
  done

done

# Remove duplicated reference sequences and rename the files
seqkit rmdup -s "$phylo_out/${results_name}_reference_sequences.fasta" \
| sed 's/ /;/' | sed 's/ /_/g' \
| awk -F';' '/^>/{print $1 ";" $(NF-1) ";" $NF; next}1' \
| sed 's/;/_/g' \
> "$phylo_out/${results_name}_reference_sequences_renamed.fasta"

cat \
"$phylo_out"/"${results_name}_high_quality_sequences.fasta" "$phylo_out/${results_name}_reference_sequences_renamed.fasta"  \
> "$phylo_out"/"${results_name}_sequences.fasta"

# Align the sequences
mafft \
  --globalpair \
  --maxiterate 1000 "$phylo_out"/"${results_name}_sequences.fasta" \
  > "$phylo_out"/"${results_name}_sequences.aln"

# Trim the sequences
trimal \
  -in "$phylo_out"/"${results_name}_sequences.aln" \
  -out "$phylo_out"/"${results_name}_sequences_trimmed.fasta" \
  -gt 0.1

# Infer the tree
iqtree2  \
  -s "$phylo_out"/"${results_name}_sequences_trimmed.fasta" \
  -bb 1000 \
  -alrt 1000 \
  -nt "$threads"

cp "$phylo_out/${results_name}_sequences_trimmed.fasta.treefile" "$results/main_outputs/${results_name}_preliminary_phylogeny.treefile"

echo ""
echo "==== PHYLOGENY COMPLETED ===="
echo "End time: $(date)"
echo "Log file saved in: $results/logs/6.phylogeny.txt"
echo ""
echo ""
