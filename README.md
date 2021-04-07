## Transcriptome annotation using nanopore sequencing

Pipeline for transcriptome annotation of genomes using nanopore sequencing technology

This project is in process, and documentation and description are still improving.

## Data preparation
To run this pipeline, following files should be prepared:
  1. Reference genome:
      * fasta: Gallus_gallus.GRCg6a.dna_sm.toplevel.fa
      * gtf: Gallus_gallus.GRCg6a.103.gtf
  2. Fastq files:
      * tissue_sampleid.fq.gz

## Configuration file
      * Chicken_Kidney_CD.fq
      * Chicken_Ovary_CC.fq
N.B. The file can not be compressed, otherwise NanoPlot would report errors

## Running the pipeline
sbatch nanopore_ann_run.sh
