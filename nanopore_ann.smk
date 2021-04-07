# Author: Dailu Guan
# Date: April 2, 2021


######################################################################################################################
##                                                                                                                  ##
##   This pipeline is used for annotating transcript isoforms using long-read sequencing of the Oxford Nanopore     ##
##   Sequencing (ONT0). The workflow is initially built by Michelle M. Halstead, who used it in the transcript      ##
##   isoform annotation of the cattle reference genome (ARS-ucd1.2).                                                ##
##                                                                                                                  ##
##   Citation: Michelle M. Halstead1, Alma Islas-Trejo1, Daniel E. Goszczynski1, Juan F. Medrano, Huaijun Zhou,     ##
##             and Pablo J. Ross. Large-scale multiplexing permits full-length transcriptome annotation of 32       ##
##             bovine tissues from a single Nanopore flow cell. (In preparation)                                    ##
##                                                                                                                  ##
######################################################################################################################

import glob
from snakemake.utils import R
import re
import socket
from time import gmtime, strftime
import pandas as pd


#### Input samples
samples_df = pd.read_csv('chicken_nanopore_samples.tsv').set_index("SampleName", drop=False)
SAMPLES = samples_df.index.values
sample = SAMPLES

#### Input chromosome
chr_df = pd.read_csv('genome_chr.txt').set_index("Chr", drop=False)
CHR = [str(i) for i in chr_df.index.values]
chr = CHR

#### input configuration ####
configfile: "config.yaml"
workdir: config['workdir']
SNAKEDIR = config['snakedir']
REF = config["ref"]
GTF = config["gtf"]


DIRS = ['01_raw_data', '02_fq_summary', '03_pychopper', '04_indiv_mapping/', '05_htseq_counts', '06_pooled_bam', '07_pooling_gff', '08_gffcompare', '09_gffread']
dir = DIRS
NANO = expand("02_fq_summary/{sample}NanoStats.txt", sample = SAMPLES)
PYCHRP = expand("03_pychopper/{sample}_pychopper_report.pdf", sample = SAMPLES)
PYCHunFQ = expand("03_pychopper/{sample}_unclassified.fq.gz", sample = SAMPLES)
PYCHrcRP = expand("03_pychopper/{sample}_rescued.fq.gz", sample = SAMPLES)
FLFP = expand("03_pychopper/{sample}_full_length.fq.gz", sample = SAMPLES)
BAM = expand("04_indiv_mapping/{sample}.bam", sample = SAMPLES)
BAI = expand("04_indiv_mapping/{sample}.bam.bai", sample = SAMPLES)
STATS = expand("04_indiv_mapping/{sample}_flagstat.txt", sample = SAMPLES)
HTCOUNT = expand("05_htseq_counts/{sample}_htseq_ENS101_counts.txt", sample = SAMPLES)
mBAM = "06_pooled_bam/allSamples_reads_aln_sorted.bam"
mBAI = "06_pooled_bam/allSamples_reads_aln_sorted.bam.bai"
mSTATS = "06_pooled_bam/allSamples_read_aln_stats.tsv"
CHRBAM = expand("07_split_bam/allSamples_reads_aln_sorted.chr{chr}.bam", chr = CHR)
CHRBAI = expand("07_split_bam/allSamples_reads_aln_sorted.chr{chr}.bam.bai", chr = CHR)
CHRGFF = expand("08_split_gff/allSamples_reads_aln_sorted.chr{chr}.gff", chr = CHR)
GFFAnn = "09_merged_gff/allSamples_reads_aln_sorted.gff"
FAnn = "11_gffread/allSamples_reads_aln_sorted_ONTann_transcriptome.fa"
MERGEFA = "11_gffread/allSamples_reads_aln_sorted_transcriptome.fa"


ruleorder: all > create_dirs > build_minimap_index > nanoplot > run_pychopper > fq_compress > indiv_mapping > run_htseq > merge_bam > split_bam > run_stringtie > merge_gffs > run_gffcompare > run_gffread

#### run the whole pipeline
rule all:
    input:
        NANO, PYCHRP, PYCHunFQ, PYCHrcRP, FLFP, BAM, BAI, STATS, HTCOUNT, mBAM, mBAI, mSTATS, CHRBAM, CHRBAI, CHRGFF, GFFAnn, FAnn, MERGEFA

#### Creating directory #####
rule create_dirs:
    output: protected(expand("{dir}", dir = DIRS))
    message: "Starting to create directory..."
    shell: "mkdir -p {dir}"


rule build_minimap_index:
    input:
        ref = REF
    output:
        index = "00_ref/minimap2_idx"
    conda:
        "Envs/minimap2.yaml"
    threads: 8
    shell:
        """
        minimap2 -t {threads} -k14 -I 1000G -d {output.index} {input.ref}
        """

rule nanoplot:
    input:
        fq = "01_raw_data/{sample}.fq"
    output:
        output = "02_fq_summary/{sample}NanoStats.txt"
    params:
        outdir = "02_fq_summary",
        prefix = "{sample}"
    conda:
        "Envs/nanoplot.yaml"
    threads: 12
    shell:
        """
        NanoPlot --fastq {input.fq} --loglength -o {params.outdir} -t {threads} -p {params.prefix} --dpi 600
        """

rule run_pychopper:
    input:
        fq = "01_raw_data/{sample}.fq"
    output:
        report = "03_pychopper/{sample}_pychopper_report.pdf",
        unclassfq = "03_pychopper/{sample}_unclassified.fq",
        rescuefq = "03_pychopper/{sample}_rescued.fq",
        flfq = "03_pychopper/{sample}_full_length.fq",
    conda:
        "Envs/pychopper.yaml"
    threads: 12
    shell:
        """
        cdna_classifier.py -t {threads} -r {output.report} -u {output.unclassfq} -w {output.rescuefq} {input} {output.flfq}
        """

rule fq_compress:
    input:
        unclassfq = "03_pychopper/{sample}_unclassified.fq",
        rescuefq = "03_pychopper/{sample}_rescued.fq",
        flfq = "03_pychopper/{sample}_full_length.fq"
    output:
        unclassfqgz = "03_pychopper/{sample}_unclassified.fq.gz",
        rescuefqgz = "03_pychopper/{sample}_rescued.fq.gz",
        flfqgz = "03_pychopper/{sample}_full_length.fq.gz"
    conda:
        "Envs/pychopper.yaml"
    threads: 8
    shell:
        """
        cat {input.unclassfq} | {SNAKEDIR}/scripts/pigz -9 -p {threads} -c > {output.unclassfqgz}
        cat {input.rescuefq} | {SNAKEDIR}/scripts/pigz -9 -p {threads} -c > {output.rescuefqgz}
        cat {input.flfq} | {SNAKEDIR}/scripts/pigz -9 -p {threads} -c > {output.flfqgz}
        """

rule indiv_mapping:
    input:
        fq = "03_pychopper/{sample}_full_length.fq.gz",
        ref = REF
    output:
        obam = "04_indiv_mapping/{sample}.bam",
        obai = "04_indiv_mapping/{sample}.bam.bai",
        state = "04_indiv_mapping/{sample}_flagstat.txt"
    conda:
        "Envs/minimap2.yaml"
    threads: 12
    shell:
        """
        minimap2 -t {threads} -ax splice -uf -k14 -G 1000000 {input.ref} {input.fq} | samtools view -q 10 -F 2304 -bS | samtools sort -@ {threads} - > {output.obam}
        samtools index -@ {threads} {output.obam}
        samtools flagstat -@ {threads} {output.obam} > {output.state}
        """

rule run_htseq:
    input:
        ibam = expand("04_indiv_mapping/{sample}.bam", sample = SAMPLES),
        ibai = expand("04_indiv_mapping/{sample}.bam.bai", sample = SAMPLES),
        gtf = GTF,
        ref = REF
    output:
        htcount = "05_htseq_counts/{sample}_htseq_ENS101_counts.txt"
    conda:
        "Envs/htseq.yaml"
    threads: 8
    shell:
        """
        htseq-count -i gene_id --format=bam --order=pos --type=exon --stranded=yes --mode=intersection-nonempty {input.ibam} {input.gtf} > {output.htcount}
        """

rule merge_bam:
    input:
        bams = expand("04_indiv_mapping/{sample}.bam", sample = SAMPLES),
        ref = REF
    output:
        obam_tmp = "06_pooled_bam/allSamples_reads_aln.bam",
        obam = "06_pooled_bam/allSamples_reads_aln_sorted.bam",
        obai = "06_pooled_bam/allSamples_reads_aln_sorted.bam.bai",
        stats = "06_pooled_bam/allSamples_read_aln_stats.tsv"
    conda:
        "Envs/samtools.yaml"
    threads: 12
    shell:
        """
        samtools merge -@ {threads} {output.obam_tmp} {input.bams}
        samtools sort -@ {threads} {output.obam_tmp} > {output.obam}
        samtools index -@ {threads} {output.obam}
        samtools flagstat -@ {threads} {output.obam} > {output.stats}
        """

rule split_bam:
    input:
        ibam = "06_pooled_bam/allSamples_reads_aln_sorted.bam",
        ibai = "06_pooled_bam/allSamples_reads_aln_sorted.bam.bai",
        ref = REF
    output:
        obam = "07_split_bam/allSamples_reads_aln_sorted.chr{chr}.bam",
        obai = "07_split_bam/allSamples_reads_aln_sorted.chr{chr}.bam.bai"
    params:
        chr = "{chr}"
    conda:
        "Envs/samtools.yaml"
    threads: 12
    shell:
        """
        samtools view -b -h -O BAM -o {output.obam} {input.ibam} {params.chr}
        samtools index -@ {threads} {output.obam}
        """

rule run_stringtie:
    input:
        bam = "07_split_bam/allSamples_reads_aln_sorted.chr{chr}.bam",
        gtf = GTF
    output:
        gff = "08_split_gff/allSamples_reads_aln_sorted.chr{chr}.gff"
    conda: "Envs/stringtie.yaml"
    threads: 12
    shell:
        """
        stringtie --rf -G {input.gtf} -l ONT -L -v -p {threads} --conservative -o {output.gff} {input.bam}
        """

rule merge_gffs:
    input:
        gffs = expand("08_split_gff/allSamples_reads_aln_sorted.chr{chr}.gff", chr = CHR)
    output:
        merged = "09_merged_gff/allSamples_reads_aln_sorted.gff"
    threads: 1
    shell: 
        """
        echo '#gff-version 2' >> {output.merged}
        echo '#pipeline-nanopore-ref-isoforms: stringtie' >> {output.merged}
        for i in {input.gffs}
        do
            grep -v '#' $i >> {output.merged}
        done
        """

rule run_gffcompare:
    input:
        ann_gff = "09_merged_gff/allSamples_reads_aln_sorted.gff",
        gtf = GTF
    output:
        cmp_dir = "10_gffcompare"
    threads: 1
    params:
        prefix = "gff_cmp"
    conda: "Envs/gffcompare.yaml"
    shell:
        """
        gffcompare -R -o {output.cmp_dir}/{params.prefix} -r {input.gtf} {input.ann_gff}
        {SNAKEDIR}/scripts/generate_tracking_summary.py --tracking 10_gffcompare/gff_cmp.tracking --output_dir {output.cmp_dir} --annotation {input.gtf}
        {SNAKEDIR}/scripts/plot_gffcmp_stats.py -r {output.cmp_dir}/gff_cmp_report.pdf -t 10_gffcompare/gff_cmp.tracking 10_gffcompare/gff_cmp.stats
        """

rule run_gffread:
    input:
        ann_gff = "09_merged_gff/allSamples_reads_aln_sorted.gff",
        ref = REF,
        gff_cmp = "10_gffcompare"
    output:
        fas = "11_gffread/allSamples_reads_aln_sorted_ONTann_transcriptome.fa",
        merged_fas = "11_gffread/allSamples_reads_aln_sorted_transcriptome.fa",
    conda: "Envs/gffread.yaml"
    threads: 1
    shell:
        """
        gffread -g {input.ref} -w {output.fas} {input.ann_gff}
        if [ -f {input.gff_cmp}/gff_cmp.annotated.gtf ]
        then
            gffread -F -g {input.ref} -w {output.merged_fas} {input.gff_cmp}/gff_cmp.annotated.gtf
        else
            touch {output.merged_fas}
        fi
        """
