#!/bin/bash

# Usage
USAGE="
Usage: ./2.trimming.sh -i <data_path> [options]

Required:
  -i  Input directory containing FASTQ files

Optional:
  -o  Output directory (default: basename of input path)
  -t  Number of threads (default: 12)

  -m  Minimum read length (default: 500)
  -M  Maximum read length (default: 2000)
  -f  Bases to trim from 5' end (default: 30)

  -w  Sliding window size for quality cutting (default: 3)
  -q  Mean quality threshold for quality cutting (default: 25)
  -p  Qualified quality phred score (default: 15)
  -Q  Mean read quality threshold (default: 20)

  -h  Show this help message
"
# Parse options
while getopts "i:o:t:m:M:f:w:q:p:Q:h" OPTION; do
  case $OPTION in
        i) data_path=$OPTARG;;
        o) results=$OPTARG;;
        t) threads=$OPTARG;;
	m) min_len=$OPTARG;;
	M) max_len=$OPTARG;;
        f) trim_front=$OPTARG;;
	w) cut_window_size=$OPTARG;;
        q) cut_mean_quality=$OPTARG;;
        p) qualified_quality_phred=$OPTARG;;
	Q) mean_qual=$OPTARG;;
        h) echo "$USAGE"; exit 0;;
        :) printf "Missing argument for option -$OPTARG\n" >&2; exit 1;;
        \?) printf "Invalid option -$OPTARG\n" >&2; exit 1;;
  esac
done

#Validate input
MISSING="is missing but required."
if [ -z "${data_path+x}" ]; then
    echo "-i $MISSING"
    echo "$USAGE"
    echo "Enjoy"
    exit 1
fi

# Set defaults
[ -z "${threads+x}" ] && threads=12
base_data_path=$(basename "$data_path")
[ -z "${results+x}" ] && results="$base_data_path"

[ -z "${min_len+x}" ] && min_len=500
[ -z "${max_len+x}" ] && max_len=2000
[ -z "${trim_front+x}" ] && trim_front=30
[ -z "${cut_window_size+x}" ] && cut_window_size=3
[ -z "${cut_mean_quality+x}" ] && cut_mean_quality=25
[ -z "${qualified_quality_phred+x}" ] && qualified_quality_phred=15
[ -z "${mean_qual+x}" ] && mean_qual=20

# Create output directory
mkdir -p "$results"

# Log
mkdir -p "$results/logs"
exec > >(tee -a "$results/logs/2.log_trimming.txt") 2>&1

#:::::::::::::::::::::::::::::::::::::::::::::::::::
# Announce configuration
echo "==== TRIMMING STARTED ===="
echo "Start time: $(date)"
echo "Input data path: $data_path"
echo "Output directory: $results/2.trimming"
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

# Trimming
output_trimmed="$results/2.trimming/trimmed"
output_fail="$results/2.trimming/fail"
output_report="$results/2.trimming/report"

mkdir -p "$output_trimmed" "$output_fail" "$output_report"

# Fastplong
echo ""
echo "=== Running fastplong trimming ==="
for file in "$data_path"/*.fastq; do
    filename=$(basename "$file")

fastplong \
        -i "$file" \
        -o "$output_trimmed/trimmed_$filename" \
        --failed_out "$output_fail/failed_trimming.$filename" \
	--length_required "$min_len" \
	--length_limit "$max_len" \
	--trim_front "$trim_front" \
	--cut_front --cut_front_window_size "$cut_window_size" --cut_front_mean_quality "$cut_mean_quality" \
	--cut_tail --cut_tail_window_size "$cut_window_size" --cut_tail_mean_quality "$cut_mean_quality" \
	--qualified_quality_phred "$qualified_quality_phred" \
	--mean_qual "$mean_qual" \
        --thread "$threads" \
        --html "$output_report/fastplong_report_$filename.html" \
        --json "$output_report/fastplong_report_$filename.json"
done
echo ""
echo "output saved in $output_trimmed"

echo ""
echo "==== TRIMMING COMPLETED ===="
echo "End time: $(date)"
echo "Log file saved in: $results/logs/2.log_trimming.txt"
echo ""
echo ""
