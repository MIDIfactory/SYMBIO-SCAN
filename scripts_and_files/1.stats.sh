#!/bin/bash

USAGE="
Usage: ./1.stats.sh -i <data_path> [options]

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
[ -z "${threads+x}" ] && threads=12
base_data_path=$(basename $data_path)
[ -z "${results+x}" ] && results="${base_data_path}_SymbiontScreener_results"
results_name=$(basename "$results")

# Create output directory
mkdir -p "$results"

# Log
mkdir -p "$results/logs"
exec > >(tee -a "$results/logs/1.log_preprocessing_QC.txt") 2>&1

#:::::::::::::::::::::::::::::::::::::::::::::::::::
# Announce configuration
echo "=== PRE PROCESSING QUALITY CHECK STARTED ===="
echo "Start time: $(date)"
echo "Input data path: $data_path"
echo "Output directory: $results/1.quality_check_pre_processing"
echo "Threads: $threads"
echo ""

# Welcoming ducks
cat <<'EOF'
 __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)< 
 \___)   \___)   \___)   \___)   \___)   \___)   \___)   \___)   
quack quack quack quack quack quack quack quack quack quack quack
EOF
echo ""

#:::::::::::::::::::::::::::::::::::::::::::::::::::

# Pre processing quality check
out_QC="$results/1.quality_check_pre_processing"
out_fastqc="$out_QC/fastqc"
out_multiqc="$out_QC/multiqc"
out_nanoplot="$out_QC/nanoplot"

mkdir -p "$out_QC" "$out_fastqc" "$out_multiqc" "$out_nanoplot"

# Run seqkit
echo ""
echo "=== Running seqkit ==="
seqkit stats -T "$data_path"/*.fastq > "$out_QC/${results_name}_stats_pre_trimming.txt"

cat "$out_QC/${results_name}_stats_pre_trimming.txt"

echo ""
echo "output saved in $out_QC/${results_name}_stats_pre_trimming.txt"

# Run fastqc
echo ""
echo "=== Running Fastqc ==="
for file in "$data_path"/*.fastq; do
    fastqc "$file" -o "$out_fastqc"
done
echo ""
echo "output saved in $out_fastqc"

# Run multiqc
echo ""
echo "=== Running Multiqc ==="
multiqc \
    "$out_fastqc"/*.zip \
    -o "$out_multiqc"
mv "$out_multiqc/multiqc_report.html" "$out_multiqc/${results_name}_multiqc_pre_report.html"
echo ""
echo "output saved in $out_multiqc/${results_name}_multiqc_pre_report.html"

# Run NanoPlot on individual samples
echo ""
echo "=== Running NanoPlot ==="
for file in "$data_path"/*.fastq; do
    filename=$(basename "$file")
    NanoPlot \
        --fastq "$data_path"/"$filename" \
        --plots dot \
        -t "$threads" \
        -o "$out_nanoplot/$filename"
done
echo -e "\noutput saved in $out_nanoplot/$filename"

# Create a folder for the more relevant results
mkdir -p "$results/main_outputs"
cp "$out_QC/${results_name}_stats_pre_trimming.txt" "$results/main_outputs/."
cp "$out_multiqc/${results_name}_multiqc_pre_report.html" "$results/main_outputs/."

echo ""
echo "==== PRE PROCESSING QUALITY CHECK ===="
echo "End time: $(date)"
echo "Log file saved in: $results/logs/1.log_preprocessing_QC.txt"
echo ""
echo ""
