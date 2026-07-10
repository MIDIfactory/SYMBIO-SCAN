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

# Create output directory
mkdir -p "$results"

# Log file
mkdir -p "$results/logs"
exec > >(tee -a "$results/logs/4.log_kraken_bracken.txt") 2>&1

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

# Run Kraken2: filtered reads (confidence 0.0)
echo ""
echo "=== Running Kraken2 on raw reads (confidence 0.0) ==="

kraken_out="$results/4.kraken/kraken_raw/kraken"
mkdir -p "$kraken_out"

for raw_read in "$data_path"/*fastq; do
    filename=$(basename "$raw_read" .fastq)
    kraken2 \
        --use-names \
        --db "$KRAKEN_DB" \
        --threads "$threads" \
        --confidence 0.0 \
        --report "$kraken_out/$filename.0.0_raw.report" \
        --classified-out "$kraken_out/$filename.classified.fq" \
        --unclassified-out "$kraken_out/$filename.unclassified.fq" \
        "$raw_read" > "$kraken_out/$filename.0.0_raw.kraken"
	
done

echo ""
echo "output saved in $kraken_out"

# Generate Krona input from Kraken reports
echo ""
echo "=== Generating Krona input ==="

krona_out="$results/4.kraken/kraken_raw/krona"
mkdir -p "$krona_out"

for report in "$kraken_out"/*.report; do
    filename=$(basename "$report" .0.0_raw.report)
    python "$SCRIPT_DIR/kreport2krona.py" -r "$report" -o "$krona_out/$filename.krona"
done

# Generate Krona HTML report
echo "=== Generating Krona HTML visualization ==="
krona_res_fix=$(echo $result | sed 's/\//_/g')
ktImportText "$krona_out"/*.krona -o "$krona_out/${results}_krona_raw.html"

echo ""
echo "output saved in $krona_out/${results}_krona_raw.html"


########################################################################
# Run Kraken2: filtered reads (confidence 0.1) and Braken at genus level

# Run Kraken2: filtered reads (confidence 0.1)
echo ""
echo "=== Running Kraken2 on filtered reads (confidence 0.1) ==="

filtered_reads="$results/2.trimming/trimmed"
kraken_out="$results/4.kraken/kraken_filt/kraken"
mkdir -p "$kraken_out"

for file in "$filtered_reads"/*fastq; do
    filename=$(basename "$file" .fastq)
    kraken2 \
        --use-names \
        --db "$KRAKEN_DB" \
        --threads "$threads" \
        --confidence 0.1 \
        --report "$kraken_out/$filename.report" \
	    --classified-out "$kraken_out/$filename.classified.fq" \
	    --unclassified-out "$kraken_out/$filename.unclassified.fq" \
        "$file" > "$kraken_out/$filename.kraken"

done

echo ""
echo "output saved in $kraken_out"


# Bracken at genus level
echo ""
echo "=== Running Braken at genus level ==="

# Define input, output and database paths
input_bracken="$kraken_out"
bracken_out="$results/4.kraken/kraken_filt/bracken"
mkdir -p "$bracken_out"

# Define taxonomic levels (available: 'D','P','C','O','F','G','S') - with SILVA only max tax level is genus
tax_levels=(G)

# Loop over Kraken2 report from filtered reads
for kraken_report in "$input_bracken"/*.report; do
    filename=$(basename "$kraken_report" .report)
    # Loop over every taxonomic level
    for level in "${tax_levels[@]}"; do
        output_file="${bracken_out}/${filename}_${level}.txt"

        # Run Bracken
        bracken \
		-d "$KRAKEN_DB" \
                -i "$kraken_report" \
                -o "$output_file" \
                -t 1 \
                -r 500 \
                -l "$level"

        # Print progress message for the current level
        echo "Processed ${kraken_report} at level ${level}, output saved to ${output_file}"
    done
done

# Print completion message
echo "All Bracken abundance estimations completed."
echo "output saved in $bracken_out"


# Generate Krona input from Kraken reports
input_report="$results/4.kraken/kraken_filt/kraken"
krona_out="$results/4.kraken/kraken_filt/krona"
mkdir -p "$krona_out"

echo ""
echo "=== Generating Krona input from Bracken reports ==="
for bracken_report in "$input_report"/*bracken*.report; do
    filename=$(basename "$bracken_report" _bracken_genuses.report | sed 's/trimmed_//g')
    python "$SCRIPT_DIR/kreport2krona.py" -r "$bracken_report" -o "$krona_out/$filename.krona"
done

# Generate Krona HTML report and open it in background
echo "=== Generating Krona HTML visualization ==="
ktImportText "$krona_out"/*.krona -o "$krona_out/${results}_krona_filt.html"

echo ""
echo "output saved in $krona_out/${results}_krona_bracken.html"



##########
# BARPLOTS

echo ""
echo "=== Producing barplots ==="

# Raw reads barplots

mkdir -p $results/4.kraken/kraken_raw/barplots

python "$SCRIPT_DIR/plot_microbiome_raw.py" \
    -i $results/4.kraken/kraken_raw/krona \
    -o $results/4.kraken/kraken_raw/barplots

echo ""
echo "output saved in $results/4.kraken/kraken_raw/barplots/stacked_barplot_top30.pdf"

# Filtered reads barplot

mkdir -p $results/4.kraken/kraken_filt/barplots

python "$SCRIPT_DIR/plot_microbiome.py" \
	 -i $results/4.kraken/kraken_filt/bracken \
	 -o $results/4.kraken/kraken_filt/barplots

echo ""
echo "output saved in $results/4.kraken/kraken_filt/barplots/stacked_barplot.pdf"

# Create a folder for the more relevant results
mkdir -p "$results/main_outputs"
cp "$results"/4.kraken/kraken_*/krona/*.html "$results/main_outputs/."
cp "$results"/4.kraken/kraken_*/barplots/*barplot* "$results/main_outputs/."
cp "$results"/4.kraken/kraken_raw/barplots/stacked_barplot_top30.pdf "$results/main_outputs/."
cp "$results"/4.kraken/kraken_filt/barplots/stacked_barplot.pdf "$results/main_outputs/."

echo ""
echo "==== KRAKEN COMPLETED ===="
echo "End time: $(date)"
echo "Log file saved in: $results/logs/4.log_kraken_bracken.txt"
echo ""
echo ""
