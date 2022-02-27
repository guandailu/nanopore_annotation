#!/bin/bash -l

# Author: Dailu Guan
# Date: April 11, 2021

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --job-name=nanopore_ann
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/group/zhougrp/dguan/nanopore_annotation/Chicken/98_logs/%x-%j.out
#SBATCH --error=/group/zhougrp/dguan/nanopore_annotation/Chicken/98_logs/%x-%j.err
#SBATCH --time=30-00:00:00
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=dguan@ucdavis.edu


#############################################
# CONDA ENV ACTIVATE
module load bio3

#############################################
# SUBMIT TASKS HERE

snakemake -j 68 \
	--cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -t {cluster.time} -n {cluster.cpus} --mem={cluster.mem} -J {cluster.name} -o {cluster.output} -e {cluster.error}" \
	--cluster-config /group/zhougrp/dguan/nanopore_annotation/Chicken/99_scripts/nanopore_ann_final.cluster.yaml \
	-s /group/zhougrp/dguan/nanopore_annotation/Chicken/99_scripts/nanopore_ann_final.2.smk \
	--configfile /group/zhougrp/dguan/nanopore_annotation/Chicken/99_scripts/config.yaml \
	--latency-wait 560 -p -k --nolock --rerun-incomplete --use-conda $@
