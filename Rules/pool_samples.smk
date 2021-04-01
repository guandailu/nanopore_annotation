import glob
from snakemake.utils import R
import re
import socket
from time import gmtime, strftime
import pandas as ph

#### Input samples
samples_df = pd.read_csv('chicken_nanopore_samples.tsv').set_index("sample", drop=False)
SAMPLES = list(samples_df['sample'])
sample = SAMPLES

#### input configuration ####
configfile: "config.yaml"
workdir: config['workdr']
REF = config["reffa"]
GTF = config["gtf"]

rule pooling_seq:
    input:
        fq = expand("01_raw_data/{sample}.fq.gz", sample = SAMPLES)
    output:
        pq = "06_pooling_seq/allSamples_reads.fq.gz",
        flq = "06_pooling_seq/allSamples_full_length_reads.fq",
        report = "06_pooling_seq/allSamples_pychopper_report.pdf",
        unclassfq = "06_pooling_seq/allSamples_unclassified.fq",
        rescuefq = "06_pooling_seq/allSamples_rescued.fq"
    params:
        pc = "True" if config["run_pychopper"] else "False",
        pc_opts = config["pychopper_opts"],
        concat = "True" if config["concatenate"] else "False",
    threads: 12
    conda: "99_scripts/Envs/pychopper.yaml"
    shell:
        """
        cat {input.fq} | gzip -c > {output.pq}
        cdna_classifier.py -t {threads} -r {output.report} -u {output.unclassfq} -w {output.rescuefq} {output.fq} {output.flq} 
        /group/zhougrp/dguan/nanopore_annotation/Chicken/99_scripts/pipeline-nanopore-ref-isoforms/scripts/generate_pychopper_stats.py --data {output.report} --output 06_pooling_seq/
        gzip {output.flq}
        gzip {output.unclassfq}
        gzip {output.rescuefq}
        """

rule pooling_mapping:
    input:
       ref = REF,
       fastq = "06_pooling_seq/allSamples_full_length_reads.fq.gz",
    output:
       bam = "07_pooling_mapping/reads_aln_sorted.bam",
       stats = "07_pooling_mapping/read_aln_stats.tsv"
    params:
    conda: '99_scripts/Envs/minimap2.yaml'
    threads: 12
    shell:"""
    minimap2 -t {threads} -ax splice -uf -k14 -G 1000000 {input.ref} {input.fastq} \
    | samtools view -q 10 -F 2304 -Sb -\
    | samtools sort -@ {threads} -o {output.bam} -;
    samtools index {output.bam}
    samtools flagstat -@ {threads} {output.bam} > {output.state}
    """

rule plot_aln_stats:
    input:
        stats = "07_pooling_mapping/read_aln_stats.tsv"
    output:
        pdf = "07_pooling_mapping/read_aln_stats.pdf"
    conda: '99_scripts/Envs/minimap2.yaml'
    threads: 12
    shell:"""
    /group/zhougrp/dguan/nanopore_annotation/Chicken/99_scripts/plot_aln_stats.py {input.stats} -r {output.pdf}
    """
