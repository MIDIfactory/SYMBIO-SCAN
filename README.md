
# SYMBIO-SCAN
SYMBIO-SCAN is a modular workflow for processing full-length 16S rRNA gene amplicons generated with Oxford Nanopore Technologies. The full bioinformatic analysis workflow was designed to process reads from each individual sample independently. The pipeline performs read quality control, trimming, taxonomic classification, consensus sequence generation and preliminary phylogenetic analysis. 

## Requirments
- **Disk space:** at least 2 GB of free space
- **CPU:** the pipeline can run with a single thread, but **at least 4 threads are recommended** for reasonable performance
- **RAM:** at least 8 GB is recommended
- **Operating system:** Linux or another Unix-like environment with Bash
- **Conda:** conda or mamba is required to install the pipeline dependencies

## Installation
The pipeline requires several tools. Most of them are compatible and can be installed together in the same Conda environment. This can be done either using the provided .yml file or by installing the packages directly.

Option 1 — using the environment file
```
git clone https://github.com/MIDIfactory/SYMBIO-SCAN.git
cd SYMBIO-SCAN

mamba env create -f environments/environment.yml
```
Option 2 — direct installation
```
mamba create -n symbio-scan -c conda-forge -c bioconda --strict-channel-priority fastqc=0.12.1 multiqc=1.31 nanoplot=1.46.1 fastplong=0.4.1 kraken2=2.1.6 krona=2.8.1 vsearch=2.30.0 minimap2=2.30 racon=1.5.0 emboss=6.6.0 trimal=1.5.0 
```
The following tools must be installed separately from the Conda environment and made available in $PATH:
```
- BLAST 2.17.0
- SeqKit 2.3.0
- Samtools 1.19.2
- MAFFT 7.505
- iqtree2 2.0.7
```
They should be installed separately from the conda environment, as they may conflict with other tools. It is recommended that you add them to your `$PATH`, otherwise the pipeline will not be able to locate them. These versions are the ones used during development and are recommended for replicating the original configuration, but **other versions should also work**.

Once the Conda environment is activated and the required tools are available in $PATH, the complete set of software used by the pipeline is:
```
fastqc	0.12.1
multiqc	1.31
NanoPlot	1.46.1
fastplong	0.4.1
kraken2	2.1.6
ktImportText	2.8.1 (part of krona)
vsearch	2.30.0
minimap2	2.30
racon	1.5.0
cons	6.6.0 (part of EMBOSS)
blastn	2.17.0
seqkit	2.3.0
samtools	1.19.2
mafft	7.505
trimal  1.5.0
iqtree2  2.0.7

```

Download the SILVA database and compile it for kraken2 and blast. NB: be sure to download the files in the "16_DB" folder.  
```
cd SYMBIO-SCAN/scripts_and_files/16S_DB

conda activate symbio-scan

kraken2-build --build --db SILVA138 --threads 6

makeblastdb -in SILVA138/data/SILVA_138.1_SSURef_NR99_tax_silva.fasta -dbtype nucl 

```

## Quick start

Remember to activate the conda environment before running the pipeline: 
```
conda activate symbio-scan
```

To verify the pipeline is working:

```
./symbio-scan.sh -h
```
the pipeline usage will show up

```
Usage: ./symbio-scan.sh -i <data_path> [options]

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
  -s  Minimum number of reads required to form a cluster in the first clusterinf (default: 10)
  -z  Minimum number of reads required to form a cluster in the second clustering (default: 15)
  -e  Allow correction outside the original consensus boundaries,
      potentially extending the consensus using OTU-supported terminal
      nucleotides
  -x  Use the reference sequence to resolve nucleotide mismatches (default: True)
  -g  Use the reference sequence to resolve gap positions (default: False)
  -G  Collapse consensus sequences assigned to the same genus by
      retaining only the largest cluster (default: True).
      Set to False to disable this filtering step.

  -h  Show this help message
```
To test the pipeline
```
./symbio-scan -i test_data -o test_data_results 
```
To analysize your reads just run 

```
./symbio-scan -i <your_data_dir> 
```
By default a folder named ``<base_data_path>_SymbiontScreener_results`` will appear in your ``working_directory`` (to specify the out_dir use the -o option).

NB: unzip your fastq files first

## Main outputs
**Quality queck report**

``stats_pre_trimming.txt``: pre-procesing seqkit quality stats report of the reads 

``multiqc_post_report.html``: pre-procesing multiqc quality check report of the reads 

``stats_post_trimming.txt``: post-procesing seqkit quality stats report of the reads 

``multiqc_pre_report.html``: post-procesing multiqc quality check report of the reads 

**Bacterial community profiling**

``test_data_results_krona_filt.html`` and ``barplot_top30_filt.pdf``: kronaplot and barplot rapresenting the bacterial taxonomic classification of the trimmed filtered reads using a confidence of 0.0

**Clustering**

``test_data_results_full_summary_list.fasta``: fasta file containing all the high-quality sequences obtained

The header is tab delimitated string having the following  fields:

``>'sample' 'cluster' 'OTUs' 'size' 'perc_id' 'tax'``


- sample: sample name
- cluster: univoque name of a sequence, and formed by the sample name followed by the cluster number for that specific sample - a cluster is basically the consensus sequence produced after that the preliminary OTUs are generated and the second clustering step occurs)
- OTUs: the OTU names of the sequences that are composing each cluster after the second clustering step
- size: number of reads composing a specific sequence
- perc_id: percentage identity of a specific sequence with its most similar reference sequence
- tax: putative taxonomy of a sequence 

es. 

```
>Bm26	Bm26_cluster_1	OTU_22;OTU_18;OTU_46;OTU_40;OTU_4;OTU_41;OTU_50;OTU_111	size=772	99.414	putative HADL01018443.69.1552 Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Anaplasmataceae;Wolbachia;Calligrapha confluens
```

``test_data_results_summary_table.tsv``: tab separated file easy to use as data sheet, where the rows are "samples names", the columns are "bacterial species", and the field cells are "fasta sequences" for that specific host-associated bacteria. 

``test_data_results_full_sus_summary_list.fasta``: fasta file containing all the sequences flagged as suspicious

``test_data_results_groups_summary.tsv``: TSV file containing detailed information about the clusters and the sequences that make them up

**Overall stats and logs**

``test_data_stats.tsv``: TSV file summarizing read/sequences counts across the entire pipeline, including input, trimmed, classified, clustered, chimera-filtered, and final reads

``LOG.txt``: comprehensive log of the modules run by the pipeline

**Preliminary phylogeny**

``test_data_results_preliminary_phylogeny.treefile``: treefile containing the preliminary taxonomy of the high-quality sequences and their reference sequences for each sample



## Pipeline description

### Quality assessment

Raw reads quality assessment is performed with FastQC (v0.12.1), multiqc (v1.31) and NanoPlot (v1.46.1). 

### Read trimming

Read preprocessing is performed with fastplong (v0.4.1) to remove adapter sequences, low-quality regions, and reads outside the expected length range. Trimming and filtering parameters are optimized for long-read amplicon sequencing, including quality-based trimming at both read ends, removal of low-quality reads, and retention of reads between 500 and 2000 bp. 
In particular the following parameters are used by default: reads shorter than 500 bp or longer than 2000 bp are discarded (--length_required 500 --length_limit 2000); the first 30 bases of each read are trimmed (--trim_front 30); sliding-window quality trimming is applied to both 5′ and 3′ ends (--cut_front --cut_front_window_size 3 --cut_front_mean_quality 25 --cut_tail --cut_tail_window_size 3 --cut_tail_mean_quality 25); bases with Phred quality scores below 15 are considered unqualified (--qualified_quality_phred 15); and reads with a mean Phred quality score below 20 are removed (--mean_qual 20).

### Post processing quality assessment

Trimmed reads quality assessment is performed as previously described. 

### Bacterial community profiling 

Taxonomic classification of trimmed reads is performed using Kraken (v2.1.6) with the SILVA reference database (SILVA 138.1). Taxonomic profiles are visualized with Krona (KronaTools v2.8.1). Only reads assigned by Kraken are retained for downstream analyses, whereas unclassified reads are discarded, as they often comprised non-specific amplification products, host-derived sequences or other sequencing artifacts. 

### Consensus generation

The Kraken classified sequences are oriented using vsearch (v2.30.0) and converted to fasta format using seqkit (v2.3.0). A preliminary chimera screening is performed with vsearch (v2.30.0) using the UCHIME_REF algorithm and the GOLD database as reference. Only sequences identified as chimeric are discarded, whereas borderline sequences are retained for downstream analyses. Using only >1000 bp reads (filtered with seqkit), sequences are clustered at 97% sequence identity using vsearch with the CLUSTER_FAST algorithm, and retaining the longest sequence of each cluster as centroid representative. Only clusters with a size of at least 10 reads are retained for downstream analyses. The classified and oriented reads, including both longer and shorter sequences than 1000 bp, are mapped using minimap2 to the centroid sequences obtained in the previous step (v2.30-r1287), followed by polishing and consensus sequence generation with racon (v1.5.0). A second chimera screening is performed on the polished consensus sequences using vsearch as previously described. Taxonomic assignment of consensus sequences is performed using BLASTn (2.17.0) against the SILVA reference database (SILVA 138.1).

After polishing and taxonomic classification, several consensus sequences displayed pairwise identities greater than 97% despite exhibiting some differences upon manual inspection. For example, chimeric fragments, short duplication or inversion, poorly aligned terminal regions containing short insertions, deletions, and mismatches were frequently observed. To reduce the impact of these artefacts on downstream analyses, a reference-guided trimming approach is implemented. To identify these regions, each consensus sequence is locally aligned using BLASTn against its closest SILVA reference sequence, as determined during the taxonomic classification step. The sequences are handled using seqkit and samtools (v1.19.2). This strategy is chosen over semiglobal alignment to better handle cases where terminal regions or internal segments exhibit low-quality or structurally complex alignments, including partial duplications, inversions, chimeric sequences and poorly aligned flanking regions. Only the best high-scoring segment pairs (HSP) per query is retained, corresponding to the most contiguous and highest-confidence alignment region, assumed to represent the true 16S rRNA core sequence. This approach effectively removes non-aligned terminal regions as well as duplicated or spurious segments resulting in HSPs with lower scores and shorter alignment lengths. To further improve sequence quality, a custom Python script (Python v3.13.7) was developed to identify and remove gap-rich and mismatch-rich regions at the sequence ends based on the reference-guided alignments. The first and last 200 bp of each sequence are scanned for stretches of consecutive gaps (>5 bp), deleting the low quality flanking regions. Then, the first and last 50 bps are screened for mismatches, trimming up to the last detected mismatch within the inspection window. This procedure precisely removed poorly aligned terminal regions while preserving as much as possible the high-confidence central portion of the alignment resulting in longer final sequences. Importantly, the selected reference is used solely as an alignment scaffold to assess the reliability of terminal regions and define trimming boundaries, and no sequence correction is performed based on the reference sequence, preserving the original nucleotide composition of each consensus sequence.

As mentioned above, after the previously discussed polishing steps, some consensus sequences had percentage identity higher than 97%. Therefore, a second clustering step is performed using vsearch as described above. Sequences belonging to the resulting clusters are gathered and grouped using custom Python scripts. For each cluster, a reference-guided consensus refinement procedure is applied. Consensus sequences assigned to the same cluster are aligned using mafft (v7.505), merged into a new consensus sequence using EMBOSS (v6.6.0.0), and taxonomically re-identified using BLASTn against the SILVA database. The newly generated consensus sequence, its closest SILVA reference sequence, and all consensus sequences contributing to the cluster are subsequently aligned. The resulting alignment is used to refine the previously obtained consensus sequence using a custom Python script. Briefly, the consensus positions differing from the SILVA reference were re-computed, using the following hierarchical criteria to distinguish genuine biological variations from putative errors. A majority rule among input sequences was always prioritized. In case of “ties”, if the same base of the SILVA reference was among the equally most abundant ones, this was preferred, except for indels. In all other cases, the first-stage read size of the input sequences was considered, and the base supported by the collectively largest read size was selected. In case of residual ties, standard IUPAC codes were used, always favoring nucleotides over gaps. Moreover, the consensus length was extended when possible. Namely, if at least one input sequence had more positions than the consensus on either end, these were added to the refined consensus, applying the criteria described above to identify the most appropriate base. This approach is improving the over quality of the sequences thanks to the reference sequence as an anchor but never introducing any changes coming from it. 

The final consensus sequences are reclassified taxonomically and screened again for chimeras as described above. Clusters with an overall reads size <15 were discarded. Because some chimeric sequences may persist through the entire analysis pipeline, the final sequences assigned to the same genus are re-evaluated. When multiple consensus sequences belong to the same genus, the sequence with the smaller size (i.e., the lower number of reads contributing to that consensus sequence) is flagged as suspicious. The high-quality sequences obtained are used for the further phylogenetic analyses. 

### Preliminary phylogeny

The obtained high-quality sequences and the reference sequences used for their taxonomical classification are then very basically aligned with mafft (v7.505), trimmed with (v1.5.rev1) and a preliminary phylogenetic tree is inferred with iqtree2 (v2.0.7). 

## Pipeline overview

<img width="1151" height="1186" alt="Supplementary Figure 1 - workflow" src="https://github.com/user-attachments/assets/db22f3b0-80a5-46d4-8faf-cb9c4cf40336" />


## Citation
If our tool is useful to you, please cite:
Croci C, Iamundo G, Baharibabadi A, Sabatelli S, Sassera D, Castelli M. SYMBIO-SCAN, a full-length 16S rRNA gene-based workflow optimized for bacterial symbiont detection and discovery from long read sequencing. 
