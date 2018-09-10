# pacbio


## Installation

Python3 needs to be in PATH.

## Example

```sh
scripts/fastq_qc_to_fasta \
	-i examples/ex2.fastq -m examples/ex1-primers.fasta -p examples/analysis/ex2


scripts/blast_reads \
    -i examples/analysis/ex2_filtered_qs85_hpf_qs30_rep2_unique_primerblasted_wsauto_unique.fasta \
    -e examples/ex2-app-exons.fasta \
    -out examples/analysis/ex2
	

```