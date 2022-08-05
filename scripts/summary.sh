#!/bin/bash

set -euo pipefail

echo -e "strain\tspecies\tmlst\tnr_contigs\tassembly_length\tN50\tcoverage_metagenomic\tcoverage_specific_metagenomic"

for file in genomes/*.fasta
do
	NAME=$(basename $file .fasta)
	KRAKEN=$(awk '$4 == "S" {print $6" "$7}' kraken_out/${NAME}_kraken2_report.txt | head -n 1)
	MLST=$(awk '{print $3}' mlst/${NAME}.tsv)
	CONTIGS=$(awk -F "\t" '$1 == "# contigs" {print $2}' quast_out/${NAME}/report.tsv)
	LENGTH=$(awk -F "\t" '$1 == "Total length" {print $2}' quast_out/${NAME}/report.tsv)
	N50=$(awk -F "\t" '$1 == "N50" {print $2}' quast_out/${NAME}/report.tsv)
	COVERAGE=$(awk -F "\t" -v N=$NAME '$1 == N {print $2}' coverage.tsv)
	COVERAGE_SPECIFIC=$(awk -F "\t" -v N=$NAME '$1 == N {print $2}' coverage_specific.tsv)
	echo -e "${NAME}\t${KRAKEN}\t${MLST}\t${CONTIGS}\t${LENGTH}\t${N50}\t${COVERAGE}\t${COVERAGE_SPECIFIC}"
done
