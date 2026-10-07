#!/bin/bash

USAGE="
Usage: ./5.amplicons_seqs.sh -i <data_path> [options]

Required:
  -i  Input directory containing FASTQ files

Optional:
  -o  Output directory (default: basename of input path)
  -l  Minimum read length to retain (default: 1000)

  -c  First OTU clustering percent identity (default: 97)
      NOTE: specify values from 0 to 100 (e.g. 97), NOT decimals (e.g. 0.97)

  -C  Second OTU clustering percent identity (default: 97)
      NOTE: specify values from 0 to 100 (e.g. 97), NOT decimals (e.g. 0.97)

  -s  Minimum number of reads required to form a cluster in the first clustering step (default: 10)

  -z  Minimum number of reads required to form a cluster in the second clustering step (default: 15)

  -e  Allow correction outside the original consensus boundaries,
      potentially extending the consensus using OTU-supported terminal
      nucleotides

  -x  Use the reference sequence to resolve nucleotide mismatches (default: True)

  -g  Use the reference sequence to resolve gap positions (default: False)

  -G  Collapse consensus sequences assigned to the same genus by
      retaining only the largest cluster (default: True)
      Set to False to disable this filtering step.


  -t  Number of threads (default: 12)

  -h  Show this help message
"

# Parse options
while getopts "i:o:l:c:C:s:z:e:x:g:G:t:h" OPTION; do
  case $OPTION in
        i) data_path=$OPTARG;;
        o) results=$OPTARG;;
        l) min_len=$OPTARG;;
        c) id_1=$OPTARG;;
        C) id_2=$OPTARG;;
        s) min_size=$OPTARG;;
        z) min_size_2=$OPTARG;;
        e) extend_ends=$OPTARG;;
        x) ref_mismatch=$OPTARG;;
        g) ref_gap=$OPTARG;;
        G) collapse_genus=$OPTARG;;
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
[ -z "${min_len+x}" ] && min_len=1000
[ -z "${id_1+x}" ] && id_1=97
perc_id_1="0.${id_1}"
[ -z "${id_2+x}" ] && id_2=97
perc_id_2="0.${id_2}"
[ -z "${min_size+x}" ] && min_size=10
[ -z "${min_size_2+x}" ] && min_size_2=15
[ -z "${extend_ends+x}" ] && extend_ends=True
[ -z "${ref_mismatch+x}" ] && ref_mismatch=True
[ -z "${ref_gap+x}" ] && ref_gap=False
[ -z "${collapse_genus+x}" ] && collapse_genus=True
[ -z "${threads+x}" ] && threads=12

# Create output directory
mkdir -p "$results"

# Log
mkdir -p "$results/logs"
exec > >(tee -a "$results/logs/5.log_amplicons_seqs.txt") 2>&1

#:::::::::::::::::::::::::::::::::::::::::::::::::::
# Announce configuration
echo "This is 5.amplicons_seqs.sh"
echo "==== HIGH QUALITY OTUs GENERATION ===="
echo "Start time: $(date)"
echo "Input data path: $data_path"
echo "Output directory: $results/5.amplicons_seqs"
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
# OTUs output dir
OTUs_out="$results/5.amplicons_seqs"
mkdir -p "$OTUs_out"

# Create the summary table output
full_groups_summary="${OTUs_out}/${results_name}_full_groups_summary.tsv"
touch "${full_groups_summary}"

full_sus_summary_list="${OTUs_out}/${results_name}_full_sus_summary_list.fasta"
touch "${full_sus_summary_list}"

full_small_size_summary_list="${OTUs_out}/${results_name}_full_small_size_summary_list.fasta"
touch "${full_sus_summary_list}"

full_summary_list="${OTUs_out}/${results_name}_full_summary_list.fasta"
touch "${full_summary_list}"

stats="${OTUs_out}/${results_name}_clustering_stats.tsv"
echo -e "sample\toriented\tINPUT_orient_nochim_reads\tinput_long_reads\t1_chim\t1_clust\t1_cons\t2_chim\t2_clust\t3_chim\tsus\tsmall_size\t2_cons" > "${stats}"

summary_table="${OTUs_out}/${results_name}_summary_table.tsv"

# Define scripts and databases
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SILVA_DB="$SCRIPT_DIR/16S_DB/SILVA138/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta"
SILVA_GOLD="$SCRIPT_DIR/16S_DB/gold.fasta"
SILVA_SUBSET="$SCRIPT_DIR/16S_DB/silva_subset_100.fasta"

### Process OTUs sample by sample
for file in "$data_path"/*.fastq; do
  sample=$(basename "$file" .fastq)
  echo ""
  echo "###########################"
  echo "Processing sample: $sample"

  # Generate sample-specific folder
  sample_out="$OTUs_out/$sample"
  mkdir -p "$sample_out"

  # 1. Orient sequences and convert fastq seqs in fasta
  echo ""
  echo "step 1. Orient sequences and convert fastq seqs in fasta"

  # Orient sequences
  classified_fq="$results/4.kraken/kraken_filt/kraken/trimmed_${sample}.classified.fq"
  oriented_out="$sample_out/1.oriented"
  mkdir -p "$oriented_out"

  vsearch \
      --orient "$classified_fq" \
      --db "$SILVA_SUBSET" \
      --fastqout "$oriented_out/${sample}_oriented.fastq" \
      --notmatched "$oriented_out/${sample}_notmatched.fastq" \
      --tabbedout "$oriented_out/${sample}_info.txt"

  # Convert fastq seqs in fasta
  echo ""
  echo "Convert fastq seqs in fasta"
  seqkit fq2fa "$oriented_out/${sample}_oriented.fastq" > "$oriented_out/${sample}_oriented.fasta"

  echo "outputs saved in $oriented_out"

  # 2. Preliminary chimera_ref detection on reads
  echo ""
  echo "step 2. Chimera detection"

  # Chimera detection
  first_chimera_out="$sample_out/2.first_ref_chimeras"
  mkdir -p "$first_chimera_out"
  vsearch \
  --uchime_ref "$oriented_out/${sample}_oriented.fasta" \
  --nonchimeras "$first_chimera_out/${sample}_oriented_nochim.fasta" \
  --chimeras "$first_chimera_out/${sample}_oriented_chimeras.fasta" \
  --borderline "$first_chimera_out/${sample}_oriented_borderline.fasta" \
  --uchimeout "$first_chimera_out/${sample}_oriented_uchime_stats.txt" \
  --db "$SILVA_GOLD" \
  --sizein

  echo "outputs saved in $first_chimera_out"

  # Merge no chimeric and borderline sequences - basically remove just chimeric sequences
  cat "$first_chimera_out/${sample}_oriented_borderline.fasta" >> "$first_chimera_out/${sample}_oriented_nochim.fasta"

  # 3. Trimming of the short reads to improve the subsequent clustering (default 1000)
  echo ""
  echo "step 3. Trimming of the short reads to improve the subsequent clustering (default 1000)"

  trim_out="$sample_out/3.trimming"
  
  mkdir -p "$trim_out"

  seqkit seq -m "${min_len}" "$first_chimera_out/${sample}_oriented_nochim.fasta" > "$trim_out/${sample}_oriented_nochim_${min_len}.fasta"

  echo "outputs saved in $trim_out"

  # 4. First clustering - generate clusters represented by the longest seq as centroids (default 97)
  echo ""
  echo "step 4. First clustering - generate clusters represented by the longest seq as centroids (default 97)"

  centroids_out="$sample_out/4.first_clustering"
  mkdir -p "$centroids_out"
  vsearch \
  --cluster_fast "$trim_out/${sample}_oriented_nochim_${min_len}.fasta" \
  --id "${perc_id_1}" \
  --relabel OTU_ \
  --sizeout \
  --threads "$threads"  \
  --centroids "$centroids_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}.fasta" \
  --uc "$centroids_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}.uc"

  echo "outputs saved in $centroids_out"

  # 5. Select only the clusters which are composed by a minimum size, which is the number of reads making a cluster (default 10)
  echo ""
  echo "step 5. Select only the clusters which are composed by a minimum size, which is the number of reads making a cluster (default 10)"

  min_size_cluster="$sample_out/5.size_filtered"
  mkdir -p "$min_size_cluster"
  vsearch \
  --sortbysize "$centroids_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}.fasta" \
  --minsize "$min_size" \
  --output "$min_size_cluster/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}.fasta"

  echo "outputs saved in $min_size_cluster"

  # 6. Mapping of the reads in the centroids (also the reads < 1000 bp)
  echo ""
  echo "step 6. Mapping of the reads in the centroids (also the reads < 1000 bp)"

  mapping_out="$sample_out/6.mapping"
  mkdir -p "$mapping_out"
  minimap2 \
  -ax map-ont \
  -K 1G \
  -t $threads \
  "$min_size_cluster/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}.fasta" \
  "$oriented_out/${sample}_oriented.fasta" \
  > "$mapping_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}.sam"

  echo "outputs saved in $mapping_out"

  # 7. Polishing and consensus generation
  echo ""
  echo "step 7. Polishing and consensus generation"

  polishing_out="$sample_out/7.polishing"
  mkdir -p "$polishing_out"
  racon \
  "$oriented_out/${sample}_oriented.fasta" \
  "$mapping_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}.sam" \
  "$min_size_cluster/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}.fasta" \
  -t $threads \
  > "$polishing_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso.fasta"

  echo "outputs saved in $polishing_out"

  # 8. Second round of chimera detection after polishing
  echo ""
  echo "step 8. Second round of chimera detection after polishing"

  second_chimera_out="$sample_out/8.second_ref_chimeras"
  mkdir -p "$second_chimera_out"
  vsearch \
  --uchime_ref "$polishing_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso.fasta" \
  --nonchimeras "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim.fasta" \
  --chimeras "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_chimeras.fasta" \
  --borderline "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_borderline.fasta" \
  --uchimeout "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_uchime_stats.txt" \
  --db "$SILVA_GOLD" \
   --sizein

  echo "outputs saved in $second_chimera_out"

  cat "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_borderline.fasta" \
  >> "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim.fasta"

  # 9. Taxonomic classification of the consensus using SILVA
  echo ""
  echo "step 9. Taxonomic classification of the consensus using SILVA"

  tax_out="$sample_out/9.taxonomy_and_blast"
  mkdir -p "$tax_out"
  blastn \
  -query "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim.fasta" \
  -db "$SILVA_DB" \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 1 \
  -num_threads $threads \
  > "$tax_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_tax.txt"

  echo "outputs saved in $tax_out"

  # 10. Reference based trimming
  echo ""
  echo "step 10. Reference based trimming"

  # Extraxt reference sequences
  ref_trim_out="$sample_out/10.reference_based_trimming"
  ref_seqs="$sample_out/10.reference_based_trimming/ref_seqs"
  mkdir -p "$ref_seqs"
      > "$ref_seqs/reference_sequences.fasta"
  cut -f2 "$tax_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_tax.txt" | \
  while IFS= read -r i; do
      samtools faidx "$SILVA_DB" "$i" \
        >> "$ref_seqs/reference_sequences.fasta"
  done

  # Retrieve the names of the reference sequences
  cut -f2,3 "$tax_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_tax.txt" | sort -u | \
  while read -r A B; do
      sed -i "s/$A/$B/g" "$ref_seqs/reference_sequences.fasta"
  done

  # Make blast database
  makeblastdb -in "$ref_seqs/reference_sequences.fasta" -dbtype nucl

  # Split multifasta containing sequences into single files to do pairwise alignment
  split_out="$ref_trim_out/splitted_seqs"
  mkdir -p "$split_out"
  seqkit split \
  -i "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim.fasta" \
  -O "$split_out"

  # Use blastn to perform a pair-wise local alignment
  blast_out="$ref_trim_out/blast_out"
  mkdir -p "$blast_out"

  for i in "$split_out"/*; do
  name=$(basename "$i")
  blastn \
  -query $i \
  -db "$ref_seqs/reference_sequences.fasta" \
  -outfmt 16 \
  -out "$blast_out/blast_out_${name}.xml";
  done

  # Trimming delle estremità
  blast_trim_out="$ref_trim_out/blast_trim_out"
  mkdir -p "$blast_trim_out"

  for xml in "$blast_out"/blast_out_*.xml; do
    name=$(basename "$xml" .xml)
    python "$SCRIPT_DIR/gaps_mismatch_trimming.py" \
        -i "$xml" \
        -o "$blast_trim_out/trimmed_${name}.fasta"
  done

  # Ricomponiamo il multifasta
  cat $blast_trim_out/trimmed_*.fasta > "$blast_trim_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_ends_trimmed.fasta"

 # Reference guided consensus generation and correction
  echo ""
  echo "step 11. Reference guided consensus generation and correction"

  
  clusters_out="$sample_out/11.second_clustering"
  cons_correction="$sample_out/12.cons_correction"
  align_out="$cons_correction/1.alignment"
  consensus_out="$cons_correction/2.consensus"
  final_otus="$cons_correction/3.final_otus"
  third_chimera_out="$sample_out/13.third_chimera_check"

  mkdir -p "$cons_correction" "$clusters_out" "$align_out" "$consensus_out" "$final_otus" "$third_chimera_out"

  # 11 Second clustering - collapse similar consensus in the same cluster
  echo ""
  echo "step 11. Second clustering - collapse similar consensus in the same cluster"

  vsearch --cluster_fast "$blast_trim_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_ends_trimmed.fasta" \
  --id "$perc_id_2" \
  --iddef 0 \
  --relabel OTU_ \
  --sizeout \
  --threads $threads \
  --centroids "$clusters_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_otus_${id_2}_consenso_nochim_otus.fasta" \
  --uc "$clusters_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_otus_${id_2}_consenso_nochim_otus.uc"

  echo "outputs saved in $clusters_out"

  # Create the file groups.tsv
  echo "running scripts_and_files/make_groups_tsv.py"

  python "$SCRIPT_DIR/make_groups_tsv.py" \
  -u "$clusters_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_otus_${id_2}_consenso_nochim_otus.uc" \
  -t "$tax_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_tax.txt" \
  -s "$sample" \
  -o "$clusters_out"

  echo ""
  echo "file ${sample}_groups.tsv saved in $clusters_out"

  # Creates multifasta containing cluster of sequences divided according to the previously defined groups
  echo ""
  echo "running scripts_and_files/sorting_cluster.py"

  python "$SCRIPT_DIR/sorting_cluster.py" \
  -i "$blast_trim_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_nochim_ends_trimmed.fasta" \
  -g "$clusters_out/${sample}_groups.tsv" \
  -s "$sample" \
  -o "$clusters_out"

  echo "files saved in $clusters_out"

  # 12. Start reference guided consensus generation and correction
  echo ""
  echo "step 12. Start reference guided consensus generation and correction"

  # Process the clusters to get trimmed consensus
  touch "$final_otus/${sample}_final_otus.fasta"

  for cluster in "$clusters_out"/*cluster*.fasta; do
    cluster_name=$(basename "$cluster" .fasta)
    nseq=$(grep -c "^>" "$cluster")
    
    # If there is only one cluster then save it in the final output
    if [ "$nseq" -eq 1 ];
     then
     sed "s/^>.*/>${cluster_name}/" "$cluster" >> "$final_otus/${sample}_final_otus.fasta"

     # Else create a consensus starting from at least two clusters
     else

     # 12.1 Align the sequences
     echo ""
     echo "step 12.1. Align the sequences"

     mafft --maxiterate 1000 --localpair "$cluster" > "$align_out/${cluster_name}.aln"
     cons -sequence "$align_out/${cluster_name}.aln" -outseq "$consensus_out/${cluster_name}_tmp_cons.fasta" -name "${cluster_name}"
     sed 's/n//g' "$consensus_out/${cluster_name}_tmp_cons.fasta" > "$consensus_out/${cluster_name}_cons.fasta"
     rm "$consensus_out/${cluster_name}_tmp_cons.fasta"

     # 12.2 consensus generation and correciton
     echo ""
     echo "step 12.2. consensus generation and correciton"
     # Taxonomic classification of the consensus
     blastn \
     -query "$consensus_out/${cluster_name}_cons.fasta" \
     -db "$ref_seqs/reference_sequences.fasta" \
     -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
     -max_target_seqs 1 \
     -num_threads $threads \
     > "$consensus_out/${cluster_name}_cons_tax.txt"

     #Take the OTUs
     cat "$align_out/${cluster_name}.aln" | sed 's/-//g' >> "${consensus_out}/sequences_${cluster_name}.fasta"

     #Take the consensus
     cat "$consensus_out/${cluster_name}_cons.fasta" >> "${consensus_out}/sequences_${cluster_name}.fasta"

     #Take the reference
     ref=$(cut -f 3 "$consensus_out/${cluster_name}_cons_tax.txt" | cut -d ' ' -f 1)
      if [ -n "$ref" ];
        then
        samtools faidx "$ref_seqs/reference_sequences.fasta" "$ref" >> "${consensus_out}/sequences_${cluster_name}.fasta"
        else
        echo "WARNING: No reference found for $cluster_name" >&2
      fi

     # Align the sequences
     mafft \
     --maxiterate 1000 \
     --localpair "${consensus_out}/sequences_${cluster_name}.fasta" \
     > "${consensus_out}/tmp_sequences_${cluster_name}.aln"
     
     seqkit sort "${consensus_out}/tmp_sequences_${cluster_name}.aln"  \
     > "${consensus_out}/sequences_${cluster_name}.aln"

     rm "${consensus_out}/tmp_sequences_${cluster_name}.aln"

     # Consensus correction
     echo "Hello I'm running scripts_and_files/cons_correction.py"
     # sequences_${cluster_name}.aln contains the reference, the consensus and the first-stage consensus making it
     python "$SCRIPT_DIR/cons_correction.py" \
     -i "${consensus_out}/sequences_${cluster_name}.aln" \
     -r "$ref" \
     -o "$consensus_out/${cluster_name}_corrected_cons.aln" \
     --extend-ends True \
     --ref-mismatch True \
     --ref-gap False

     cat "${consensus_out}/sequences_${cluster_name}.aln" "$consensus_out/${cluster_name}_corrected_cons.aln" \
     > "${consensus_out}/sequences_${cluster_name}_corrected.aln"

     # remove gaps
     sed 's/-//g' "$consensus_out/${cluster_name}_corrected_cons.aln" > "$consensus_out/${cluster_name}_corrected_cons.fasta"
     cat "$consensus_out/${cluster_name}_corrected_cons.fasta" >> "$final_otus/${sample}_final_otus.fasta"
     rm "$consensus_out/${cluster_name}_corrected_cons.aln"

    fi
  done

  # 12.3 Taxonomic classification of the final OTUs using SILVA
  echo ""
  echo "step 12.3 Taxonomic classification of the polished consensus sequences"

  sed -i 's/_corrected//g' "$final_otus/${sample}_final_otus.fasta"
  blastn \
  -query "$final_otus/${sample}_final_otus.fasta" \
  -db "$ref_seqs/reference_sequences.fasta" \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -max_target_seqs 1 \
  -num_threads $threads \
  > "$final_otus/${sample}_final_otus_tax.txt"

  echo "outputs saved in $align_out, $consensus_out and $final_otus"

  # 13. Third round of chimera detection after polishing
  echo ""
  echo "step 13. Third round of chimera detection after polishing"

  vsearch \
  --uchime_ref "$final_otus/${sample}_final_otus.fasta" \
  --nonchimeras "$third_chimera_out/${sample}_final_otus_nochim.fasta" \
  --chimeras "$third_chimera_out/${sample}_final_otus_chimeras.fasta" \
  --borderline "$third_chimera_out/${sample}_final_otus_borderline.fasta" \
  --uchimeout "$third_chimera_out/${sample}_final_otus_uchime_stats.txt" \
  --db "$SILVA_GOLD" \
  --sizein

  # Keep borderline chimeric sequences
  cat "$third_chimera_out/${sample}_final_otus_borderline.fasta" >> "$third_chimera_out/${sample}_final_otus_nochim.fasta"

  # 14. Gather everything in the final results folder and rename the files
  echo ""
  echo "step 14. Gather everything in the final results folder and rename the files"

  final_seqs="$sample_out/14.final_sequences"
  mkdir "$final_seqs"

  # Rename the final OTUs file
  echo "running scripts_and_files/final_otus_header_rename.py"

  python "$SCRIPT_DIR/final_otus_header_rename.py" \
  -i "$third_chimera_out/${sample}_final_otus_nochim.fasta" \
  -t "$final_otus/${sample}_final_otus_tax.txt" \
  -g "$clusters_out/${sample}_groups.tsv" \
  -s "$sample" \
  -o "$final_seqs" \
  -G "$collapse_genus" \
  -z "$min_size_2"

  echo "Renamed final OTUs written to $final_seqs/${sample}_high_quality_seqs.fasta"
  echo "Small-size sequences written to $final_seqs/${sample}_small_size_sequences.fasta"
  echo "Sus sequences written to $final_seqs/${sample}_sus_sequences.fasta"
  
  # Join the results from different samples in one table easy to read in excel
  echo "Concatenating"
  cat "$clusters_out/${sample}_groups.tsv" >> "${full_groups_summary}"
  cat "$final_seqs/${sample}_high_quality_seqs.fasta" >> "${full_summary_list}"
  cat "$final_seqs/${sample}_small_size_sequences.fasta" >> "${full_small_size_summary_list}"
  cat "$final_seqs/${sample}_sus_sequences.fasta" >> "${full_sus_summary_list}"

  # Calculate the stats 
  ORIENT=$(grep -c ">" "$oriented_out/${sample}_oriented.fasta")
  INPUT_orient_nochim_reads=$(grep -c ">" "$first_chimera_out/${sample}_oriented_nochim.fasta")
  INPUT_long_reads=$(grep -c ">" "$trim_out/${sample}_oriented_nochim_${min_len}.fasta")
  CHIM_1=$(grep -c ">" "$first_chimera_out/${sample}_oriented_chimeras.fasta")
  CLUST_1=$(grep -c ">" "$centroids_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}.fasta")
  CONS_1=$(grep -c ">" "$polishing_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso.fasta")
  CHIM_2=$(grep -c ">" "$second_chimera_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_consenso_chimeras.fasta")
  CLUST_2=$(grep -c ">" "$clusters_out/${sample}_oriented_nochim_${min_len}_otus_${id_1}_min${min_size}_otus_${id_2}_consenso_nochim_otus.fasta")
  CHIM_3=$(grep -c ">" "$third_chimera_out/${sample}_final_otus_chimeras.fasta")
  SUS=$(grep -c ">" "$final_seqs/${sample}_sus_sequences.fasta")
  SMALL_SIZE=$(grep -c ">" "$final_seqs/${sample}_small_size_sequences.fasta")
  CONS_2=$(grep -c ">" "$final_seqs/${sample}_high_quality_seqs.fasta")

  echo -e "${sample}\t${ORIENT}\t${INPUT_orient_nochim_reads}\t${INPUT_long_reads}\t${CHIM_1}\t${CLUST_1}\t${CONS_1}\t${CHIM_2}\t${CLUST_2}\t${CHIM_3}\t${SUS}\t${SMALL_SIZE}\t${CONS_2}" >> "${stats}"


done

echo "running scripts_and_files/summary_table.py"
python "$SCRIPT_DIR/summary_table.py" \
-i "${full_summary_list}" \
-o "${summary_table}"

# Create a folder for the more relevant results
echo "copying results"
mkdir -p "$results/main_outputs"
cp "${full_groups_summary}" "$results/main_outputs/."
cp "${full_sus_summary_list}" "$results/main_outputs/."
cp "${full_summary_list}" "$results/main_outputs/."
cp "${summary_table}" "$results/main_outputs/."

echo ""
echo "==== OTUs GENERATION COMPLETED ===="
echo "End time: $(date)"
echo "Log file saved in: $results/logs/5.log_amplicons_seqs.txt"
echo ""
echo ""
