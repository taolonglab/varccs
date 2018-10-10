# VarCCS: Variant Calling Tool for Circular Consensus Sequencing

To do:

Removed 1 or 2 hardcoded places in blast_reads referening to first and last exons.
Automatically figure out first and the last from the exon reference
FASTA file.


## Installation

Python environment present (either 2.7 or 3+).

Required packages for R: data.table, Biostrings, stringr, optparse, ggplot2, circlize, RColorBrewer.
Required packages for Python: .

## Example

In this example we will analyze DNA sequences obtained from prefrontal cortex neuron DNA, using PCR to amplify any sequences that share the beginning and the end with the coding sequence of an Amyloid Precursor Protein (APP).

These amplified sequences were sequenced using PacBio SMRT platform and the raw BAM were converted to Circular Consensus (CCS) FASTQ files using a CCS2 algorithm in the SMRTLink version 4.0.

To start, we need several files:

1. Input CCS FASTQ sequences; here `examples/ex1.fastq`.
2. Reference exon sequences; here `examples/ex1-primers.fasta` [(tutorial how to
prepare this file)](prepare_reference_exons.md)
3. Prepared table with Familial AD mutations `examples/ex1-app-fads.txt`, if we
want to analyze them (not required).

```sh
# Quality Control 
scripts/fastq_qc_to_fasta \
	-i examples/ex1.fastq -m examples/ex1-primers.fasta -p examples/analysis/ex1
```
```
FASTQ Quality Control to FASTA

Jobs to do:
* average quality score filter: ON (min qscore: 85)
* homopolymer run filter: ON (qscore threshold: 30, num. repeats: 2+)
* PCR primer alignment filter: ON (BLAST word size: auto)

Executing jobs:
* average quality filter... OK.
* homopolymer run filter... OK.
* PCR primer alignment filter...
  ? word size automatically set to 3/4 of the shortest primer length.
  - BLAST word size: 13, gap open: 0, gap extend: 2, max offset: 150 nt
  - total reads in the read file: 257. file format: FASTA
  - processed   256 out of   257 reads (skipped:   45,   45 missing > 1 primer,  113 rev. complemented,  0 wrong dir).
  - total run time: 0 m 10 s
  OK.
* finding unique reads (adding up counts)...  OK.
  -> examples/analysis/ex1_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_wsauto_unique.fasta
* deleting temporary files...... OK.

Done.
```

In this pipeline, the output files are indicated with `->` symbol in the script outputs.

```sh
# Blast reads vs reference exons
scripts/blast_reads \
    -i examples/analysis/ex1_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_wsauto_unique.fasta \
    -e examples/ex1-app-exons.fasta \
    --prefix examples/analysis/ex1

```
```
BLAST Quality Controlled reads to Reference Exons

* calculating read lengths... OK.
* blast reads vs ref. exons
  - word size: 25, gap open: 0, gap extend: 2
  - processed   175 out of   175 reads
  - total run time: 0 m 8 s
* calculating reference exon lengths... OK.
* BLAST results analysis
  - loading and processing data...OK.
  - generating various plots...OK.
  - calculating exon statistics...OK.
  - writing exon combinations... OK.
    -> examples/analysis/ex1_exon_combinations_summary.txt
  - writing final table... OK.
* Fixing low quality homopolymer runs
  - loading BLAST file... OK.
  - loading FASTA file... OK.
  - correcting reads... OK.
  - saving output FASTA... OK.
    -> examples/analysis/ex1_hmpfix.fasta
  - saving output homopolymer fixed table... OK.
    -> examples/analysis/ex1_exon_table_hmpfix.txt
* deleting temporary files..... OK.

Done.
```

```sh
# Create a table with exon-exon join information
Rscript scripts/exon_joins.R \
    -i examples/analysis/ex1_exon_table_hmpfix.txt \
    -o examples/analysis/ex1_exon-exon_joins.txt
```

```
Exon-exon join analysis
* Loading and processing data... OK.
  - using homopolymer corrected sequences
* Finding and organizing all exon-exon joins... OK.
* saving exon-joins table... OK.
  -> examples/analysis/ex1_exon-exon_joins.txt
```

```sh
# Run SNV/indel analysis
scripts/snvs_indels_analysis \
    -i examples/analysis/ex1_exon_table_hmpfix.txt \
    -e examples/ex1-app-exons.fasta \
    -f examples/ex1-app-fads.txt \
    -os examples/analysis/ex1_snvs_indels.txt \
    -of examples/analysis/ex1_fads.txt
```

```
SNVs and indels analysis

* Generating exon lengths... OK.
* Initializing R libraries... OK.
* Loading input files..examples/ex1-app-fads.txt. OK.
* Searching for SNVs and indels... OK.
* Saving SNV/indel table... OK.
  -> examples/analysis/ex1_snvs_indels.txt
* Analyzing data for Familial AD mutations... OK.
* Saving identified known variants... OK.
  -> examples/analysis/ex1_fads.txt
* Deleting temporary files... OK.

Done.
```

```
# Circo plot 1: SNVs, indels and intra-exon joins
scripts/snvs_indels_plot \
    -i examples/analysis/ex1_exon_table_hmpfix.txt \
    -e examples/ex1-app-exons.fasta \
    -s examples/analysis/ex1_snvs_indels.txt \
    -o examples/analysis/ex1_snvs_indels_circo.pdf
```

```
SNVs and indels plots

* Generating exon lengths... OK.
* Loading R libraries... OK.
* Loading files.. OK.
* Combining SNVs and indels into a single table... OK.
* Finding all intra-exon joins... OK.
* Generating circo plot... OK.
* Deleting temporary files... OK.
```

```
# Run reading frame analysis
scripts/reading_frame_analysis \
	-i examples/analysis/ex1_exon_table_hmpfix.txt \
	-e examples/ex1-app-exons.fasta \
	-o examples/analysis/ex1_reading_frame_contigs.txt
```

```
Reading frame analysis

* Loading fixed BLAST results... OK.
* Loading exon ORF information...OK.
* Filtering out possible PCR chimeras...OK.
* Calculating positions for reverse-complemented reads...OK.
* Generating table of reading frame contiguous segments...OK.
* Savint output table... OK.
  -> examples/analysis/ex1_reading_frame_contigs.txt
* Deleting temporary files... OK.
```

```
# Plot circo reading frames and exon composition
Rscript scripts/plot_circo_reading_frames.R \
  --input-contigs examples/analysis/ex1_reading_frame_contigs.txt \
  --input-exons examples/analysis/ex1_exon_table_hmpfix.txt \
  --input-fasta examples/analysis/ex1_hmpfix.fasta \
  --output-circoplot examples/anaylsis/ex1_reading_frames_circo.pdf

```