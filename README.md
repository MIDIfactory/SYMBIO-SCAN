# SYMBIONT-SCREENER
A bioinformatic approach to obtain high quality sequences by full-length 16S amplicons 

## 1. Overview

GuanoPipeline is a modular workflow for processing full-length 16S rRNA gene amplicons generated with Oxford Nanopore Technologies. The pipeline was designed to obtain a general snapshot of the bacterial community associated with a host of interest, as well as high-quality sequences of those symbionts. The pipeline performs read quality control, trimming, taxonomic classification, consensus sequence generation and preliminary phylogenetic analysis. 

<img width="793" height="275" alt="Screenshot from 2026-07-10 15-38-46" src="https://github.com/user-attachments/assets/01e163a0-0d4f-4888-951e-c15638b81f8b" />

## 2. Installation
THIS SECTION IS STILL ON GOING - DO NOT FOLLOW IT
```
https://github.com/MIDIfactory/SYMBIONT-SCREENER.git
cd GuanoPipeline
conda env create -f environments/environment.yml
conda activate guanopipeline
```

## 3. Quick start

To analysize your reads just run 

```
./guanopipeline -i <your_data_dir> 
```

By default a folder named ``<base_data_path>_SymbiontScreener_results`` will appear in your ``working_directory`` (to specify the out_dir use the -o option).


Using the command 

```
./guano_pipeline.sh -h
```
the pipeline usage will show up

```
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
  -I  OTU clustering percent identity (default: 97)
      NOTE: specify values from 0 to 100 (e.g. 97), NOT decimals (e.g. 0.97)
  -s  Minimum number of reads required to form a cluster (default: 10)
  -e  Allow correction outside the original consensus boundaries,
      potentially extending the consensus using OTU-supported terminal
      nucleotides
  -x  Use the reference sequence to resolve nucleotide mismatches (default: True)
  -g  Use the reference sequence to resolve gap positions (default: False)
  -c  Collapse consensus sequences assigned to the same genus by
      retaining only the largest cluster (default: True).
      Set to False to disable this filtering step.

  -h  Show this help message
```

## 4. Main outputs
**Quality queck report**

``stats_pre_trimming.txt``: pre-procesing seqkit quality stats report of the reads 

``multiqc_post_report.html``: pre-procesing multiqc quality check report of the reads 

``stats_post_trimming.txt``: post-procesing seqkit quality stats report of the reads 

``multiqc_pre_report.html``: post-procesing multiqc quality check report of the reads 

**Bacterial community profiling**

``test_data_results_krona_raw.html`` and ``stacked_barplot_top30.pdf``: kronaplot and barplot rapresenting the bacterial taxonomic classification of the untrimmed raw reads using a confidence of 0.0 (aka giving a taxonomic classification whenever possible to get as broad as possible picture of the reads composition for each sample)

``test_data_results_krona_filt.html`` and ``stacked_barplot.png``: kronaplot and barplot rapresenting the bacterial taxonomic classification of the trimmed filtered reads using a confidence of 0.1 (aka giving a more conservative taxonomic classification of the reads composition for each sample)

**High-quality sequences**

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


``test_data_results_full_sus_summary_list.fasta``: fasta file containing all the suspicious sequences obtained

``test_data_results_summary_table.tsv``: tab separated file easy to use as data sheet, where the rows are "samples names", the columns are "bacterial species", and the field cells are "fasta sequences" for that specific host-associated bacteria. 

**Rough preliminary phylogeny**

``test_data_results_preliminary_phylogeny.treefile``: treefile containing the preliminary taxonomy of the high-quality sequences and their reference sequences for each sample



## 5. Pipeline description

### 5.1 Quality assessment

Raw reads quality assessment is performed with FastQC (v0.12.1), multiqc (v1.31) and NanoPlot (v1.46.1). 

### 5.2 Read trimming

Read preprocessing is performed with fastplong (v0.4.1) to remove adapter sequences, low-quality regions, and reads outside the expected length range. Trimming and filtering parameters are optimized for long-read amplicon sequencing, including quality-based trimming at both read ends, removal of low-quality reads, and retention of reads between 500 and 2000 bp. 
In particular the following parameters are used by default: reads shorter than 500 bp or longer than 2000 bp are discarded (--length_required 500 --length_limit 2000); the first 30 bases of each read are trimmed (--trim_front 30); sliding-window quality trimming is applied to both 5′ and 3′ ends (--cut_front --cut_front_window_size 3 --cut_front_mean_quality 25 --cut_tail --cut_tail_window_size 3 --cut_tail_mean_quality 25); bases with Phred quality scores below 15 are considered unqualified (--qualified_quality_phred 15); and reads with a mean Phred quality score below 20 are removed (--mean_qual 20).

### 5.3 Post processing quality assessment

Trimmed reads quality assessment is performed as previously described. 

### 5.4 Bacterial community profiling 

Taxonomic classification of trimmed reads is performed using Kraken (v2.1.6) with the SILVA reference database (SILVA 138.1). Genus-level abundance estimates are subsequently refined using Bracken (v3.0.1). Taxonomic profiles are visualized with Krona (KronaTools v2.8.1). Only reads assigned by Kraken are retained for downstream analyses, whereas unclassified reads are discarded, as they often comprised non-specific amplification products, host-derived sequences or other sequencing artifacts. 

### 5.5 Consensus generation

The Kraken classified sequences are oriented using vsearch (v2.30.0) and converted to fasta format using seqkit (v2.3.0). A preliminary chimera screening is performed with vsearch (v2.30.0) using the UCHIME_REF algorithm and the GOLD database as reference. Only sequences identified as chimeric are discarded, whereas borderline sequences are retained for downstream analyses. Length filtering is subsequently applied using seqkit to retain sequences at least 1000 bp long. 

Sequences are clustered at 97% sequence identity using vsearch with the CLUSTER_FAST algorithm, retaining the longest sequence of each cluster as centroid representative. Only clusters with a size of at least 10 reads are retained for downstream analyses. The classified and oriented reads, including both longer and shorter sequences than 1000 bp, are mapped using minimap2 to the centroid sequences obtained in the previous step (v2.30-r1287), followed by polishing and consensus sequence generation with racon (v1.5.0). A second chimera screening is performed on the polished consensus sequences using vsearch as previously described. Taxonomic assignment of consensus sequences is performed using BLASTn (2.17.0) against the SILVA reference database (SILVA 138.1).

After polishing and taxonomic classification, several consensus sequences displayed pairwise identities greater than 97% despite exhibiting some differences upon manual inspection. For example, chimeric fragments, short duplication or inversion, poorly aligned terminal regions containing short insertions, deletions, and mismatches were frequently observed. To reduce the impact of these artefacts on downstream analyses, a reference-guided trimming approach is implemented. To identify these regions, each consensus sequence is locally aligned using BLASTn against its closest SILVA reference sequence, as determined during the taxonomic classification step. The sequences are handled using seqkit and samtools (v1.19.2). This strategy is chosen over semiglobal alignment to better handle cases where terminal regions or internal segments exhibit low-quality or structurally complex alignments, including partial duplications, inversions, chimeric sequences and poorly aligned flanking regions. Only the best high-scoring segment pairs (HSP) per query is retained, corresponding to the most contiguous and highest-confidence alignment region, assumed to represent the true 16S rRNA core sequence. This approach effectively removes non-aligned terminal regions as well as duplicated or spurious segments resulting in HSPs with lower scores and shorter alignment lengths. To further improve sequence quality, a custom Python script (Python v3.13.7) was developed to identify and remove gap-rich and mismatch-rich regions at the sequence ends based on the reference-guided alignments. The first and last 200 bp of each sequence are scanned for stretches of consecutive gaps (>5 bp), deleting the low quality flanking regions. Then, the first and last 50 bps are screened for mismatches, trimming up to the last detected mismatch within the inspection window. This procedure precisely removed poorly aligned terminal regions while preserving as much as possible the high-confidence central portion of the alignment resulting in longer final sequences. Importantly, the selected reference is used solely as an alignment scaffold to assess the reliability of terminal regions and define trimming boundaries, and no sequence correction is performed based on the reference sequence, preserving the original nucleotide composition of each consensus sequence.

As mentioned above, after the previously discussed polishing steps, some consensus sequences had percentage identity higher than 97%. Therefore, a second clustering step is performed using vsearch as described above. Sequences belonging to the resulting clusters are retrieved and grouped using custom Python scripts. For each cluster, a reference-guided consensus refinement procedure is applied. Consensus sequences assigned to the same cluster are aligned using mafft (v7.505) (Katoh and Standley 2013), merged into a new consensus sequence using EMBOSS (v6.6.0.0), and taxonomically re-identified using BLASTn against the SILVA database. The newly generated consensus sequence, its closest SILVA reference sequence, and all consensus sequences contributing to the cluster are subsequently aligned. The resulting alignment is used to refine the previously obtained consensus sequence using a custom Python script. Briefly, only positions where the newly generated consensus differed from the SILVA reference are considered potentially erroneous and subjected to further inspection. For each such position, the nucleotide composition of all contributing sequences is evaluated. The correction procedure followed a hierarchical decision framework. First, the presence of a majority nucleotide among sequences is assessed based on the number of sequences supporting each nucleotide state, irrespective of their abundance. When a nucleotide is supported by more than half of the sequences, the consensus position is corrected accordingly. In the absence of a majority, the reference nucleotide is adopted whenever at least one sequence supported the same state. In this second step gap positions are treated independently, and reference gaps are ignored by default as insertions and deletions could represent genuine biological variation. If neither condition is met, sequence abundance is used as a confidence metric, and the nucleotide supported by the most abundant sequence is selected. When equally abundant sequences supported different nucleotides, standard IUPAC ambiguity codes are used to represent the uncertainty, whereas ties involving a nucleotide and a gap are resolved in favour of the nucleotide. This approach is improving the over quality of the sequences thanks to the reference sequence as an anchor but never introducing any changes coming from it.

The final consensus sequences are reclassified taxonomically and screened again for chimeras as described above. Because some chimeric sequences may persist through the entire analysis pipeline, the final sequences assigned to the same genus are re-evaluated. When multiple consensus sequences belong to the same genus, the sequence with the smaller size (i.e., the lower number of reads contributing to that consensus sequence) is flagged as a potential chimera. This approach is based on the assumption that the co-occurrence of multiple bacterial taxa from the same genus within a sample is less likely than the presence of a chimeric artifact. The high-quality sequences obtained are used for the further phylogenetic analyses. 

### 5.6 Preliminary phylogeny

The obtained high-quality sequences and the reference sequences used for their taxonomical classification are then very basically aligned with mafft (v7.505), trimmed with (v1.5.rev1) and a preliminary phylogenetic tree is inferred with iqtree2 (v2.0.7). 
