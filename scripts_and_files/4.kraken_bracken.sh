#!/bin/bash

# USAGE message
USAGE="
Usage: ./4.kraken_bracken.sh -i <data_path> [options]

Required:
  -i  Input directory containing FASTQ files

Optional:
  -o  Output directory (default: basename of input path)
  -t  Number of threads (default: 12)

  -h  Show this help message
"

# Parse options
while getopts "i:o:t:h" OPTION; do
  case $OPTION in
    i) data_path=$OPTARG ;;
    o) results=$OPTARG ;;
    t) threads=$OPTARG ;;
    h) echo "$USAGE"; exit 0 ;;
    :) printf "Missing argument for option -$OPTARG\n" >&2; echo "$USAGE"; exit 1 ;;
    \?) printf "Invalid option: -$OPTARG\n" >&2; echo "$USAGE"; exit 1 ;;
  esac
done

# Validate required input
MISSING="is missing but required."
if [ -z "${data_path+x}" ]; then
    echo "-i $MISSING"
    echo "$USAGE"
    echo "Enjoy"
    exit 1
fi

# Set defaults
[ -z "${threads+x}" ] && threads=12
base_data_path=$(basename $data_path)
[ -z "${results+x}" ] && results="${base_data_path}_SymbiontScreener_results"
results_name=$(basename "$results")

# Create output directory
mkdir -p "$results"

# Log file
mkdir -p "$results/logs"
exec > >(tee -a "$results/logs/4.log_kraken.txt") 2>&1

#:::::::::::::::::::::::::::::::::::::::::::::::::::
# Announce configuration
echo ""
echo "==== KRAKEN TAXONOMIC ASSIGNMENT STARTED ===="
echo "Start time: $(date)"
echo "Input path: $data_path"
echo "Output directory: $results"
echo ""

# Welcoming ducks
cat <<'EOF'
 __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<
 \___)   \___)   \___)   \___)   \___)   \___)   \___)   \___)
quack quack quack quack quack quack quack quack quack quack quack
EOF
echo ""
#:::::::::::::::::::::::::::::::::::::::::::::::::::

# Define scripts and Kraken DB path
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
KRAKEN_DB="$SCRIPT_DIR/16S_DB/kraken_db"

#########################################
# Run Kraken2: raw reads (confidence 0.0)
echo ""
echo "=== Running Kraken2 on raw reads (confidence 0.0) ==="

kraken_out_raw="$results/4.kraken/kraken_raw/kraken"
mkdir -p "$kraken_out_raw"

kraken_raw_stats="${kraken_out_raw}/${results_name}_kraken_raw_stats.tsv"
echo -e "sample\tunclassified_raw\tclassified_raw\tP_raw\tC_raw\tO_raw\tF_raw\tG_raw" > "${kraken_raw_stats}"

for raw_read in "$data_path"/*fastq; do
    filename=$(basename "$raw_read" .fastq)

   kraken2 \
        --use-names \
        --db "$KRAKEN_DB" \
        --threads "$threads" \
        --confidence 0.0 \
        --report "$kraken_out_raw/$filename.0.0_raw.report" \
        --classified-out "$kraken_out_raw/$filename.classified.fq" \
        --unclassified-out "$kraken_out_raw/$filename.unclassified.fq" \
        "$raw_read" > "$kraken_out_raw/$filename.0.0_raw.kraken"
	
    UNC_raw=$(grep -c ">" "$kraken_out_raw/$filename.unclassified.fq")
    CLAS_raw=$(grep -c ">" "$kraken_out_raw/$filename.classified.fq")
    P_raw=$(grep -c -w P "$kraken_out_raw/$filename.0.0_raw.report")
    C_raw=$(grep -c -w C "$kraken_out_raw/$filename.0.0_raw.report")
    O_raw=$(grep -c -w O "$kraken_out_raw/$filename.0.0_raw.report")
    F_raw=$(grep -c -w F "$kraken_out_raw/$filename.0.0_raw.report")
    G_raw=$(grep -c -w G "$kraken_out_raw/$filename.0.0_raw.report")
    echo -e "${filename}\t${UNC_raw}\t${CLAS_raw}\t${P_raw}\t${C_raw}\t${O_raw}\t${F_raw}\t${G_raw}" >> "${kraken_raw_stats}"

done

echo ""
echo "output saved in $kraken_out_raw"

# Generate Krona input from Kraken reports
echo ""
echo "=== Generating Krona input ==="

krona_out_raw="$results/4.kraken/kraken_raw/krona"
mkdir -p "$krona_out_raw"

for report in "$kraken_out_raw"/*.report; do
    filename=$(basename "$report" .0.0_raw.report)
    python "$SCRIPT_DIR/kreport2krona.py" -r "$report" -o "$krona_out_raw/$filename.krona"
done

# Generate Krona HTML report
echo "=== Generating Krona HTML visualization ==="
krona_res_fix=$(echo $result | sed 's/\//_/g')
ktImportText "$krona_out_raw"/*.krona -o "$krona_out_raw/${results_name}_krona_raw.html"

echo ""
echo "output saved in $krona_out_raw/${results_name}_krona_raw.html"


########################################################################
# Run Kraken2: filtered reads (confidence 0.0) 

echo ""
echo "=== Running Kraken2 on filtered reads (confidence 0.0) ==="

filtered_reads="$results/2.trimming/trimmed"
kraken_out_filt="$results/4.kraken/kraken_filt/kraken"
mkdir -p "$kraken_out_filt"

kraken_filt_stats="${kraken_out_filt}/${results_name}_kraken_filt_stats.tsv"
echo -e "sample\tunclassified_filt\tclassified_filt\tP_filt\tC_filt\tO_filt\tF_filt\tG_filt" > "${kraken_filt_stats}"

for file in "$filtered_reads"/*fastq; do
    filename=$(basename "$file" .fastq)
    kraken2 \
        --use-names \
        --db "$KRAKEN_DB" \
        --threads "$threads" \
        --confidence 0.0 \
        --report "$kraken_out_filt/$filename.report" \
	    --classified-out "$kraken_out_filt/$filename.classified.fq" \
	    --unclassified-out "$kraken_out_filt/$filename.unclassified.fq" \
        "$file" > "$kraken_out_filt/$filename.kraken"

    UNC_filt=$(grep -c ">" "$kraken_out_filt/$filename.unclassified.fq")
    CLAS_filt=$(grep -c ">" "$kraken_out_filt/$filename.classified.fq")
    P_filt=$(grep -c -w P "$kraken_out_filt/$filename.report")
    C_filt=$(grep -c -w C "$kraken_out_filt/$filename.report")
    O_filt=$(grep -c -w O "$kraken_out_filt/$filename.report")
    F_filt=$(grep -c -w F "$kraken_out_filt/$filename.report")
    G_filt=$(grep -c -w G "$kraken_out_filt/$filename.report")
    echo -e "${filename}\t${UNC_filt}\t${CLAS_filt}\t${P_filt}\t${C_filt}\t${O_filt}\t${F_filt}\t${G_filt}" >> "${kraken_filt_stats}"
done

echo ""
echo "output saved in $kraken_out_filt"

# Generate Krona input from Kraken reports
echo ""
echo "=== Generating Krona input ==="

krona_out_filt="$results/4.kraken/kraken_filt/krona"
mkdir -p "$krona_out_filt"

for report in "$kraken_out_filt"/*.report; do
    filename=$(basename "$report" .report)
    python "$SCRIPT_DIR/kreport2krona.py" -r "$report" -o "$krona_out_filt/$filename.krona"
done

# Generate Krona HTML report
echo "=== Generating Krona HTML visualization ==="
krona_res_fix=$(echo $result | sed 's/\//_/g')
ktImportText "$krona_out_filt"/*.krona -o "$krona_out_filt/${results_name}_krona_filt.html"

echo ""
echo "output saved in $krona_out_filt/${results_name}_krona_filt.html"

##########
# BARPLOTS

echo ""
echo "=== Producing barplots ==="

# Raw reads barplots

mkdir -p $results/4.kraken/kraken_raw/barplots

python "$SCRIPT_DIR/plot_microbiome_raw.py" \
    -i $results/4.kraken/kraken_raw/krona \
    -o $results/4.kraken/kraken_raw/barplots

mv "$results/4.kraken/kraken_raw/barplots/stacked_barplot_top30.pdf" "$results/4.kraken/kraken_raw/barplots/barplot_top30_raw.pdf"

echo ""
echo "output saved in $results/4.kraken/kraken_raw/barplots/barplot_raw.pdf"

# Filtered reads barplot

mkdir -p $results/4.kraken/kraken_filt/barplots

python "$SCRIPT_DIR/plot_microbiome_raw.py" \
    -i $results/4.kraken/kraken_filt/krona \
    -o $results/4.kraken/kraken_filt/barplots

mv "$results/4.kraken/kraken_filt/barplots/stacked_barplot_top30.pdf" "$results/4.kraken/kraken_filt/barplots/barplot_top30_filt.pdf"

echo ""
echo "output saved in $results/4.kraken/kraken_filt/barplots/barplot_top30_filt.pdf"

# Create a stats file

# Stats
KRAKEN_stats="$results/4.kraken/${results_name}_kraken_stats.tsv"
paste <(cat ${kraken_raw_stats}) <(cat ${kraken_filt_stats}) \
| cut -f 1,2,3,4,5,6,7,8,10,11,12,13,14,15,16 \
> "${KRAKEN_stats}"

# Create a folder for the more relevant results
mkdir -p "$results/main_outputs"
cp "$results"/4.kraken/kraken_filt/krona/*.html "$results/main_outputs/."
cp "$results/4.kraken/kraken_filt/barplots/barplot_top30_filt.pdf" "$results/main_outputs/."

echo ""
echo "==== KRAKEN COMPLETED ===="
echo "End time: $(date)"
echo "Log file saved in: $results/logs/4.log_kraken.txt"
echo ""
echo ""
