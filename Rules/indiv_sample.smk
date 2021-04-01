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


rule nanoplot:
    input: 
        fq = "01_raw_data/{sample}.fq.gz"
    output: 
        sum = "02_fq_summary/{sample}"
    params:
        outdir = "02_fq_summary"
    conda:
        '99_scripts/Envs/nanoplot.yaml'
    threads: 12
    shell:
        """
        Nanoplot --fastq {input.fq} --loglength -o {params.outdir} -t {threads} -p {output.sum} -f pdf
        """

rule pychopper:
    input: 
        fq = "01_raw_data/{sample}.fq.gz"
    output:
        report = "03_pychopper/{sample}_pychopper_report.pdf",
        unclassfq = "03_pychopper/{sample}_unclassified.fq",
        rescuefq = "03_pychopper/{sample}_rescued.fq",
        flfq = "03_pychopper/{sample}_full_length.fq"
    conda:
        '99_scripts/Envs/pychopper.yaml'
    threads: 12
    shell:
        """
        cdna_classifier.py -t {threads} -r {output.report} -u {output.unclassfq} -w {output.rescuefq} {input} {output.flfq}
        gzip {output.unclassfq}
        gzip {output.rescuefq}
        gzip {output.flfq}
        """


rule build_minimap_index:
    input:
        ref = REF
    output:
        index = "00_ref/minimap2_idx"
    conda: 
        '99_scripts/Envs/minimap2.yaml'
    threads: 12
    shell:"""
        minimap2 -t {threads} -k14 -I 1000G -d {output.index} {input.ref}
    """

rule mapping_by_sample:
    input: 
        fq = "03_pychopper/{sample}_full_length.fq.gz",
        ref = REF
    output:
        obam = "04_mapping/{sample}.bam",
        obai = "04_mapping/{sample}.bam.bai",
        state = "04_mapping/{sample}_flagstat.txt"
    conda:
        '99_scripts/Envs/minimap2.yaml'
    threads: 12
    shell:
        """
        minimap2 -t {threads} -ax splice -uf -k14 -G 1000000 {input.ref} {input.fq} | samtools view -q 10 -F 2304 -bS | samtools sort -@ {threads} - > {output.obam}
        samtools index -@ {threads} {output.obam}
        samtools flagstat -@ {threads} {output.obam} > {output.state}
        """

rule htseq:
    input: 
        ibam = expand("04_mapping/{sample}.bam", sample = SAMPLES),
        ibai = expand("04_mapping/{sample}.bam.bai", sample = SAMPLES),
        gtf = GTF,
        ref = REF
    output:
        htcount = "05_htseq_counts/{sample}_htseq_ENS101_counts.txt"
    conda:
        '99_scripts/Envs/htseq.yaml'
    threads: 1
    shell:
        """
        htseq-count -i gene_id --format=bam --order=pos --type=exon --stranded=yes --mode=intersection-nonempty {input.ibam} {input.gtf} > {output.htcount}
        """

