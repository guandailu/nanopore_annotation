#!/bin/bash -l

module load bio3

snakemake -j 24 --cluster-config /group/zhougrp/dguan/nanopore_annotation/Chicken/99_scripts/cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -n {cluster.cpus} --mem={cluster.mem} -J {cluster.name} -o {cluster.output} -e {cluster.error}" -s nanopore_ann.smk --configfile config.yaml --latency-wait 560 -p -k --use-conda $@

