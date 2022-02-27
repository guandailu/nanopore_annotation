## Full-length transcript prediction using nanopore sequencing

Pipeline for transcript prediction of the chicken genome using nanopore sequencing technology

## Data preparation
To run this pipeline, following files should be prepared:
  1. Reference genome:
      * fasta: Gallus_gallus.GRCg6a.dna_sm.toplevel.fa
      * gtf: Gallus_gallus.GRCg6a.102.gtf
  2. Fastq files:
      * tissue_sampleid.fq.gz

## Configuration file
      * Chicken_Kidney_CD.fq
      * Chicken_Ovary_CC.fq
N.B. The file can not be compressed, otherwise NanoPlot would report errors

## Running the pipeline
sbatch nanopore_ann_run.sh

N.B. The pipeline is currently configured on UC DAVIS FARM server, which will be updated in the future.
