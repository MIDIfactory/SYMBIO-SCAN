#!/bin/bash

USAGE="
Usage: ./guano_pipeline.sh -i <data_path> [options]

Required:
  -i  Input directory containing FASTQ files

General options:
  -o  Output directory (default: basename of input path)
  -t  Number of threads used by all modules (default: 12)

  -r  Comma-separated list of modules to run (default: 1,2,3,4,5,6)
      Available modules:
        1  Reads quality check
        2  Reads trimming
        3  Post-trimming quality check
        4  Taxonomic classification (Kraken2/Bracken)
        5  Consensus sequence generation
        6  Preliminary phylogeny

Module 2: Read trimming
  -m  Minimum read length after trimming (default: 500)
  -M  Maximum read length (default: 2000)
  -f  Bases to trim from the 5' end (default: 30)
  -w  Sliding window size for quality trimming (default: 3)
  -q  Mean quality threshold for quality trimming (default: 25)
  -p  Qualified quality Phred score (default: 15)
  -Q  Mean read quality threshold (default: 20)

Module 5: Consensus sequence generation
  -l  Minimum read length to retain for consensus generation (default: 1000)
  -c  First OTU clustering percent identity (default: 97)
      NOTE: specify values from 0 to 100 (e.g. 97), NOT decimals (e.g. 0.97)
  -C  Second OTU clustering percent identity (default: 97)
      NOTE: specify values from 0 to 100 (e.g. 97), NOT decimals (e.g. 0.97)
  -s  Minimum number of reads required to form a cluster (default: 10)
  -e  Allow correction outside the original consensus boundaries,
      potentially extending the consensus using OTU-supported terminal
      nucleotides
  -x  Use the reference sequence to resolve nucleotide mismatches (default: True)
  -g  Use the reference sequence to resolve gap positions (default: False)
  -G  Collapse consensus sequences assigned to the same genus by
      retaining only the largest cluster (default: True).
      Set to False to disable this filtering step.

  -h  Show this help message
"

# Parse options
while getopts "i:o:t:r:m:M:f:w:q:p:Q:l:c:C:s:e:x:g:G:h" OPTION; do
  case $OPTION in
        i) data_path=$OPTARG;;
        o) results=$OPTARG;;
        t) threads=$OPTARG;;
        r) modules=$OPTARG;;

        # Module 2: trimming
        m) min_len_trim=$OPTARG;;
        M) max_len=$OPTARG;;
        f) trim_front=$OPTARG;;
        w) cut_window_size=$OPTARG;;
        q) cut_mean_quality=$OPTARG;;
        p) qualified_quality_phred=$OPTARG;;
        Q) mean_qual=$OPTARG;;

        # Module 5: consensus generation
        l) min_len=$OPTARG;;
        c) id_1=$OPTARG;;
        C) id_2=$OPTARG;;
        s) min_size=$OPTARG;;
        e) extend_ends=$OPTARG;;
        x) ref_mismatch=$OPTARG;;
        g) ref_gap=$OPTARG;;
        G) collapse_genus=$OPTARG;;

        h) echo "$USAGE"; exit 0;;
        :) printf "Missing argument for option -$OPTARG\n" >&2; exit 1;;
        \?) printf "Invalid option -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required."

if [ -z "${data_path+x}" ]; then
    echo "-i $MISSING"
    echo "$USAGE"
    exit 1
fi

# Set defaults
[ -z "${threads+x}" ] && threads=12
base_data_path=$(basename "$data_path")
[ -z "${results+x}" ] && results="${base_data_path}_SymbiontScreener_results"
results_name=$(basename "$results")
[ -z "${modules+x}" ] && modules="1,2,3,4,5,6"

# Module 2 defaults
[ -z "${min_len_trim+x}" ] && min_len_trim=500
[ -z "${max_len+x}" ] && max_len=2000
[ -z "${trim_front+x}" ] && trim_front=30
[ -z "${cut_window_size+x}" ] && cut_window_size=3
[ -z "${cut_mean_quality+x}" ] && cut_mean_quality=25
[ -z "${qualified_quality_phred+x}" ] && qualified_quality_phred=15
[ -z "${mean_qual+x}" ] && mean_qual=20

# Module 5 defaults
[ -z "${min_len+x}" ] && min_len=1000
[ -z "${id_1+x}" ] && id_1=97
[ -z "${id_2+x}" ] && id_2=97
[ -z "${min_size+x}" ] && min_size=10
[ -z "${extend_ends+x}" ] && extend_ends=True
[ -z "${ref_mismatch+x}" ] && ref_mismatch=True
[ -z "${ref_gap+x}" ] && ref_gap=False
[ -z "${collapse_genus+x}" ] && collapse_genus=True

run_module() {
    [[ ",$modules," == *",$1,"* ]]
}

#::::::::::::::::::::::::::::::::::::::::::::::
# WELCOME
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
echo -e "\nHello world, welcome to the GuanoPipeline \nIt’s $current_date_time \nSomewhere Sun is shining, and rain is raining"

cat <<'EOF'

 __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<
 \___)   \___)   \___)   \___)   \___)   \___)   \___)   \___)
quack quack quack quack quack quack quack quack quack quack quack

EOF
echo ""

#:::::::::::::::::::::::::::::::::::::::::::::::::::

# Define module scripts position
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Run the pipeline
if run_module 1; then
bash "$SCRIPT_DIR/scripts_and_files/1.stats.sh" \
    -i "$data_path" \
    -o "$results" \
    -t "$threads"
fi

if run_module 2; then
bash "$SCRIPT_DIR/scripts_and_files/2.trimming.sh" \
    -i "$data_path" \
    -o "$results" \
    -t "$threads" \
    -m "$min_len_trim" \
    -M "$max_len" \
    -f "$trim_front" \
    -w "$cut_window_size" \
    -q "$cut_mean_quality" \
    -p "$qualified_quality_phred" \
    -Q "$mean_qual"
fi

if run_module 3; then
bash "$SCRIPT_DIR/scripts_and_files/3.quality_check_post_processing.sh" \
    -i "$data_path" \
    -o "$results" \
    -t "$threads"
fi

if run_module 4; then
bash "$SCRIPT_DIR/scripts_and_files/4.kraken_bracken.sh" \
    -i "$data_path" \
    -o "$results" \
    -t "$threads"
fi

if run_module 5; then
bash "$SCRIPT_DIR/scripts_and_files/5.amplicons_seqs.sh" \
    -i "$data_path" \
    -o "$results" \
    -l "$min_len" \
    -c "$id_1" \
    -C "$id_2" \
    -s "$min_size" \
    -e "$extend_ends" \
    -x "$ref_mismatch" \
    -g "$ref_gap" \
    -G "$collapse_genus" \
    -t "$threads"
fi

if run_module 6; then
bash "$SCRIPT_DIR/scripts_and_files/6.phylogeny.sh" \
    -i "$data_path" \
    -o "$results" \
    -t "$threads"
fi

# Generate stats
paste \
   <(sed -E 's/ +/\t/g' "$results/1.quality_check_pre_processing/${results_name}_stats_pre_trimming.txt" | cut -f 1 | rev | cut -f 1 -d '/' | rev | sed 's/.fastq//g') \
   <(sed -E 's/ +/\t/g' "$results/1.quality_check_pre_processing/${results_name}_stats_pre_trimming.txt" | cut -f 4-8) \
   <(sed -E 's/ +/\t/g' "$results/3.quality_check_post_processing/${results_name}_stats_post_trimming.txt" | cut -f 4-8) \
   <(cut -f 2- "$results/4.kraken/${results_name}_kraken_stats.tsv") \
   <(cut -f 2- "$results/5.amplicons_seqs/${results_name}_clustering_stats.tsv") \
   > "$results/main_outputs/${results_name}_stats.tsv"

# Generate logs
log_dir="$results/logs"
cat "$log_dir"/*log*txt > "$log_dir/LOG.txt"
echo ""
echo "Log file saved in: $results/logs/LOG.txt"

# Greetings
current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
echo ""
echo "Bye Bye from guano_pipeline"
echo "It’s $current_date_time"
echo "Somewhere Sun is still shining, and rain is still raining"
echo ""

cat <<'EOF'
 __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)<  __(.)< 
 \___)   \___)   \___)   \___)   \___)   \___)   \___)   \___)   
quack quack quack quack quack quack quack quack quack quack quack
EOF
