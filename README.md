# pacbio

To do:

Removed 1 or 2 hardcoded places in blast_reads referening to exons 1 and 18.
Automatically figure out first and the last from the exon reference
FASTA file.


## Installation

Python environment present (either 2.7 or 3+).

Required packages for R: data.table, Biostrings, stringr, optparse, ggplot2, circlize, RColorBrewer.
Required packages for Python: .

## Example

To start we need several files:

1. Input CCS FASTQ sequences; here `examples/ex1.fastq`.
2. Reference exon sequences; here `examples/ex1-primers.fasta` [(tutorial how to
prepare this file)](prepare_reference_exons.md)
3. Prepared table with Familial AD mutations `examples/ex1-app-fads.txt`, if we
want to analyze them (not required).

```sh
scripts/fastq_qc_to_fasta \
	-i examples/ex1.fastq -m examples/ex1-primers.fasta -p examples/analysis/ex1
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
  - total run time: 0 m 9 s
  OK.
* finding unique reads (adding up counts)...  OK.
  -> examples/analysis/ex1_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_wsauto_unique.fasta
* deleting temporary files...... OK.

Done.

```sh
# Blast reads vs reference exons
scripts/blast_reads \
    -i examples/analysis/ex1_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_wsauto_unique.fasta \
    -e examples/ex1-app-exons.fasta \
    --prefix examples/analysis/ex1

# Create a table with exon-exon join information
# I don't think this is even used in the followup analysis.
Rscript scripts/exon_joins.R \
    -i examples/analysis/ex1_exon_table_hmpfix.txt \
    -o examples/analysis/ex1_exon-exon_joins.txt

scripts/snvs_indels_analysis \
    -i examples/analysis/ex1_exon_table_hmpfix.txt \
    -e examples/ex1-app-exons.fasta \
    -os examples/analysis/ex1_snvs_indels.txt \

# Circo plot 1: SNVs, indels and intra-exon joins
scripts/snvs_indels_plot \
    -i examples/analysis/ex1_exon_table_hmpfix.txt \
    -e examples/ex1-app-exons.fasta \
    -s examples/analysis/ex1_snvs_indels.txt \
    -f examples/analysis/ex1_fads.txt \
    -o examples/analysis/ex1_snvs_indels_circo.pdf

# Run reading frame analysis
scripts/reading_frame_analysis \
	-i examples/analysis/ex1_exon_table_hmpfix.txt \
	-e examples/ex1-app-exons.fasta \
	-o examples/analysis/ex1_reading_frame_contigs.txt
```