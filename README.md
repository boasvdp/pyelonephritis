# Pyelonephritis

## Pipeline

First, the pipeline trims isolate sequencing data and assembles this using the Shovill wrapper around SPAdes. Kraken2 is used to detect species in the isolate sequencing data and Quast is used to evaluate the draft genome assembly metrics.

The draft genome is annotated in various ways:

1. Virulence genes are detected using ABRicate with the VFDB database
2. Resistance genes are detected using AMRfinderplus
3. MLST is assessed using mlst
4. General annotation using Prokka

Next, the trimmed metagenomic reads are mapped onto the draft genome using Snippy. Reads that could not be mapped are written to a separate file and checked for species using Kraken2.

Results are collected and summarised using basic bash/awk scripts.

An overview of the pipeline:
![rulegraph](https://github.com/boasvdp/pyelonephritis/blob/master/rulegraph.svg "Rulegraph")
