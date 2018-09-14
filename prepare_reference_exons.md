# Preparing reference exon sequences: APP example

## Downloading sequences from UCSC Genome Browser

We start by downloading the FASTA file with the full list of all exons and introns
for the 'app' gene:

1. Go to UCSC Genome Browser: http://genome.ucsc.edu/cgi-bin/hgGateway
2. Choose GRCh38 genome assembly, then search for APP. The automatic selection
should automatically show "Current position: chr21:25,880,550-26,170,654.".
3. Click GO, then another screen hide everything except for NCBI RefSeq on
display. We see several genes with slightly different number of Exons. For this
we will take the one with 18 exons (click on the second one fromt he top).
4. This one has the RefSeq: NM_000484.3, which will be shown on the next page.
Click on 'Genomic Sequence from assembly' link somewhere in the middle.
5. Select 5' UTR exons, CDS Exons, 3' UTR exons, introns, one FASTA per region,
and UNcheck split. Keep all exons as uppercase, everything else as lowercase.
Click 'get DNA'. Copy paste into text editor and save as
[app_exons_introns.fa](fa/app_exons_introns.fa).
 - Alternatively, we only look at CDS (so uncheck 5' and 3' UTR exons). This file
   ends with _cds_only.
6. Convert this FASTA to single line:
```sh
awk '/^>/ { printf("\n%s\n", $0);next } { printf("%s", $0) }' app_exons_introns_cds_only.fa > app_exons_introns_cds_only_singleline.fa
```

Check the exon count:
```sh
awk '/^[ACGT]/ { print $0 }' app_exons_introns.fa | wc -l
      18
```

## Extract exons and introns into separate files

This gets us the exons:
```sh
grep -B1 "^[ACGT]" app_exons_introns_cds_only_singleline.fa | grep -v -- "^--$" > app_exons_cds.fa
```

This will also get rid of the annoying group separator '--' with B1 option.
Let's also cleanup the metadata for this fasta file:

```sh
awk 'BEGIN { i=1 } /^>/ { split($0, meta, " "); sub(">hg38_ncbiRefSeqCurated_", "", meta[1]); printf(">exon_%02d %s %s\n", i, meta[1], meta[2]); i++; next } { print $0 }' app_exons_cds.fa > app_exons_cds_clean.fa
```

Now let's do introns:

```sh
grep -B1 "^[acgt]" app_exons_introns.fa | grep -v -- "^--$" > app_introns.fa
```

Clean-up:
```sh
awk 'BEGIN { i=1 } /^>/ { split($0, meta, " "); sub(">hg38_ncbiRefSeqCurated_", "", meta[1]); printf(">intron_%02d %s %s\n", i, meta[1], meta[2]); i++; next } { print $0 }' app_introns.fa > app_introns_clean.fa
```

## Generate reverse complement sequences

Run the script:
```sh
scripts/fasta_add_revcomp app_introns_clean.fa ex1-app-exons.fasta
```

We can now remove any intermediate files manually if we wish.