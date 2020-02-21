#!/bin/bash

# Goal of this script is to clean the fastq file names in a reproducible manner

for R1 in isolates/raw_reads/*R1*fq.gz
do
	R1CLEAN=$(echo ${R1%-1solate*}_R1.fastq.gz | tr -d '-')
	echo "Renaming ${R1} to ${R1CLEAN}"
	mv "$R1" "$R1CLEAN"
done

for R2 in isolates/raw_reads/*R2*fq.gz
do
	R2CLEAN=$(echo ${R2%-1solate*}_R2.fastq.gz | tr -d '-')
	echo "Renaming ${R2} to ${R2CLEAN}"
	mv "$R2" "$R2CLEAN"
done

for meta in metagenomics/*fq
do
	echo "gzipping $meta"
	gzip $meta
done
