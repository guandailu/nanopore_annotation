# Date: April 2, 2021


######################################################################################################################
##                                                                                                                  ##
##   This pipeline is used for annotating transcript isoforms using long-read sequencing of the Oxford Nanopore     ##
##   Sequencing (ONT). The workflow is initially built by Michelle M. Halstead, who used it in the transcript      ##
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

#### input configuration ####
configfile: "config.yaml"
workdir: config["workdir"]
SNAKEDIR = config["snakedir"]
REF = config["ref"]
GTF = config["gtf"]
name_prefix = config["name_prefix"]

#### Input samples
sample_file = config["SampleList"]
samples_df = pd.read_csv(sample_file).set_index("SampleName", drop=False)
SAMPLES = samples_df.index.values
sample = SAMPLES

#### Input chromosome
genome_chr = config["ChrList"]
chr_df = pd.read_csv(genome_chr).set_index("Chr", drop=False)
CHR = [str(i) for i in chr_df.index.values]
chr = CHR


NANO = expand("02_fq_summary/{sample}NanoStats.txt", sample = SAMPLES)
PYCHRP = expand("03_pychopper/{sample}_pychopper_report.pdf", sample = SAMPLES)
PYCHunFQ = expand("03_pychopper/{sample}_unclassified.fq.gz", sample = SAMPLES)
PYCHrcRP = expand("03_pychopper/{sample}_rescued.fq.gz", sample = SAMPLES)
FLFP = expand("03_pychopper/{sample}_full_length.fq.gz", sample = SAMPLES)
SAM = expand("04_indiv_mapping/{sample}.sam", sample = SAMPLES)
BAM = expand("04_indiv_mapping/{sample}.bam", sample = SAMPLES)
BAI = expand("04_indiv_mapping/{sample}.bam.bai", sample = SAMPLES)
STATS = expand("04_indiv_mapping/{sample}_flagstat.txt", sample = SAMPLES)
HTCOUNT = expand("05_htseq_counts/{sample}_htseq_ENS102_counts.txt", sample = SAMPLES)
FQ = "06_mergered_fq/allSamples.fq"
aFLFQ = "06_mergered_fq/allSamples_full_length.fq"
mBAM = "07_mergered_alignment/allSamples.bam"
mBAM = "07_mergered_alignment/allSamples_sort.bam",
mSTATES = "07_mergered_alignment/allSamples_sort.tsv"
GTFAnn = "08_stringtie/allSamples.gtf"
TRACK = "09_gffcompare/gff_cmp.tracking"
COMGFF = "09_gffcompare/gff_cmp.annotated.gtf"
MERGEFA = "10_gffread/gff_cmp.annotated_transcriptome.fa"
AL = "11_suppa/allSamples_ann_AL_strict.gtf"
aBAM = expand("12_ann_mapping/{sample}.bam", sample = SAMPLES)
TSE = expand("13_tse/{sample}.tsv", sample = SAMPLES)
CPPRED = "14_cppred/gff_cmp.annotated_transcriptome.cppred.result"


ruleorder: all > build_minimap_index > nanoplot > run_pychopper > fq_compress > indiv_mapping > run_htseq > merge_fq > preprocess_pychopper > generate_fq_stats > map_reads > clean_bam > plot_aln_stats > run_stringtie > run_gffcompare > run_gffread > run_suppa > mapping_ann > run_cppred

#### run the whole pipeline
rule all:
    input:
        NANO, PYCHRP, PYCHunFQ, PYCHrcRP, FLFP, SAM, BAM, BAI, STATS, HTCOUNT, mBAM, mSTATES, GTFAnn, TRACK, COMGFF, MERGEFA, AL, aBAM, TSE, CPPRED, aFLFQ

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
        flfq = "03_pychopper/{sample}_full_length.fq"
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
        cat {input.unclassfq} | pigz -9 -p {threads} -c > {output.unclassfqgz}
        cat {input.rescuefq} | pigz -9 -p {threads} -c > {output.rescuefqgz}
        cat {input.flfq} | pigz -9 -p {threads} -c > {output.flfqgz}
        """

rule indiv_mapping:
    input:
        fq = "03_pychopper/{sample}_full_length.fq.gz",
        ref = REF
    output:
        osam = "04_indiv_mapping/{sample}.sam",
        obam = "04_indiv_mapping/{sample}.bam",
        obai = "04_indiv_mapping/{sample}.bam.bai",
        stats = "04_indiv_mapping/{sample}_flagstat.txt"
    conda:
        "Envs/minimap2.yaml"
    threads: 12
    shell:
        """
        minimap2 -t {threads} -ax splice -uf -k14 -G 1000000 {input.ref} {input.fq} > {output.osam}
        samtools view -q 10 -F 2304 -bS {output.osam} | samtools sort -@ {threads} - > {output.obam}
        samtools index -@ {threads} {output.obam}
        samtools flagstat -@ {threads} {output.obam} > {output.stats}
        """

rule run_htseq:
    input:
        ibam = "04_indiv_mapping/{sample}.bam",
        ibai = "04_indiv_mapping/{sample}.bam.bai",
        gtf = GTF,
        ref = REF
    output:
        htcount = "05_htseq_counts/{sample}_htseq_ENS102_counts.txt"
    conda:
        "Envs/htseq.yaml"
    threads: 8
    shell:
        """
        htseq-count -i gene_id --format=bam --order=pos --type=exon --stranded=yes --mode=intersection-nonempty {input.ibam} {input.gtf} > {output.htcount}
        """


rule merge_fq:
    input:
        fqs = expand("01_raw_data/{sample}.fq", sample = SAMPLES),
        ref = REF
    output:
        fq = "06_mergered_fq/allSamples.fq",
    threads: 1
    shell:
        """
        cat {input.fqs} > {output.fq}
        """

rule preprocess_pychopper:
    input:
        fq = "06_mergered_fq/allSamples.fq"
    output:
        report = "06_mergered_fq/allSamples_pychopper_report.pdf",
        unclassfq = "06_mergered_fq/allSamples_unclassified.fq",
        rescuefq = "06_mergered_fq/allSamples_rescued.fq",
        flfq = "06_mergered_fq/allSamples_full_length.fq"
    conda:
        "Envs/pychopper.yaml"
    threads: 12
    shell:
        """
        cdna_classifier.py -t {threads} -r {output.report} -u {output.unclassfq} -w {output.rescuefq} {input} {output.flfq}
        """

rule generate_fq_stats:
    input:
        fastq = "06_mergered_fq/allSamples_full_length.fq"
    output:
        lengths = "06_mergered_fq/lengths.txt",
        base_qual = "06_mergered_fq/base_qual.txt",
        read_qual = "06_mergered_fq/read_qual.txt"
    threads: 1
    shell:
        """
        {SNAKEDIR}/scripts/run_fastq_qc.py --fastq {input.fastq} --output "06_mergered_fq"
        """


rule map_reads:
    input:
        index = "00_ref/minimap2_idx",
        fastq = "06_mergered_fq/allSamples_full_length.fq"
    output:
        sam = "07_mergered_alignment/allSamples.sam",
        bam = "07_mergered_alignment/allSamples.bam",
    conda: "Envs/minimap2.yaml"
    threads: 12
    shell:
        """
        minimap2 -t {threads} -ax splice -uf {input.index} {input.fastq} > {output.sam};
        samtools view -q 10 -F 2304 -Sb {output.sam} > {output.bam}
        """

rule clean_bam:
    input:
        bam = "07_mergered_alignment/allSamples.bam"
    output:
        sort = "07_mergered_alignment/allSamples_sort.bam",
        stats = "07_mergered_alignment/allSamples_sort.tsv"
    conda:
        "Envs/minimap2.yaml"
    threads: 12
    shell:
        """
        seqkit bam -j {threads} -x -T '{{Yaml: "99_scripts/AlnContext.yaml"}}' {input.bam} | samtools sort -@ {threads} -o {output.sort} -;
        samtools index {output.sort}
        ((seqkit bam -s -j {threads} {output.sort} 2>&1)  | tee {output.stats} ) || true
        if [[ -s 07_mergered_alignment/internal_priming_fail.tsv ]];
            then
                tail -n +2 07_mergered_alignment/internal_priming_fail.tsv | awk '{{print ">" $1 "\\n" $4 }}' - > 07_mergered_alignment/context_internal_priming_fail_start.fasta
                tail -n +2 07_mergered_alignment/internal_priming_fail.tsv | awk '{{print ">" $1 "\\n" $6 }}' - > 07_mergered_alignment/context_internal_priming_fail_end.fasta
        fi
        """

rule plot_aln_stats:
    input:
        stats = "07_mergered_alignment/allSamples_sort.tsv",
    output:
        pdf = "07_mergered_alignment/allSamples_sort.pdf"
    threads: 1
    shell:
        """
        {SNAKEDIR}/scripts/plot_aln_stats.py {input.stats} -r {output.pdf}
        """

rule run_stringtie:
    input:
        bam = "07_mergered_alignment/allSamples_sort.bam",
        gtf = GTF
    output:
        gff = "08_stringtie/allSamples.gtf"
    params:
        prefix = name_prefix
    conda: "Envs/stringtie.yaml"
    threads: 12
    shell:
        """
        stringtie --rf -G {input.gtf} -l {params.prefix} -L -v -p {threads} --conservative -o {output.gff} {input.bam}
        """

rule run_gffcompare:
    input:
        ann_gff = "08_stringtie/allSamples.gtf",
        gtf = GTF
    output:
        track = "09_gffcompare/gff_cmp.tracking",
        cmp_gff = "09_gffcompare/gff_cmp.annotated.gtf"
    threads: 1
    params:
        outdir = "09_gffcompare",
        prefix = "09_gffcompare/gff_cmp"
    conda: "Envs/gffcompare.yaml"
    shell:
        """
        gffcompare -R -o {params.prefix} -r {input.gtf} {input.ann_gff}
        {SNAKEDIR}/scripts/generate_tracking_summary.py --tracking 09_gffcompare/gff_cmp.tracking --output_dir {params.outdir} --annotation {input.gtf}
        {SNAKEDIR}/scripts/plot_gffcmp_stats.py -r {params.outdir}/gff_cmp_report.pdf -t 09_gffcompare/gff_cmp.tracking 09_gffcompare/gff_cmp.stats
        """

rule run_gffread:
    input:
        ref = REF,
        cmp_gff = "09_gffcompare/gff_cmp.annotated.gtf"
    output:
        merged_fa = "10_gffread/gff_cmp.annotated_transcriptome.fa",
    conda: "Envs/gffread.yaml"
    threads: 1
    shell:
        """
        gffread -F -g {input.ref} -w {output.merged_fa} {input.cmp_gff}
        """

rule run_suppa:
    input:
        ann_gff = "09_gffcompare/gff_cmp.annotated.gtf",
        gtf = GTF
    output:
        AL = "11_suppa/allSamples_ann_AL_strict.gtf"
    conda: "Envs/suppa.yaml",
    params:
        suppa = "11_suppa",
        ref_prefix = "Gallus_gallus_GRCg6a",
        ann_prefix = "allSamples_ann"
    threads: 1
    shell:
        """
        suppa.py generateEvents -i {input.gtf} -o {params.suppa}/{params.ref_prefix} -f ioe -e SE SS MX RI FL
        suppa.py generateEvents -i {input.ann_gff} -o {params.suppa}/{params.ann_prefix} -f ioe -e SE SS MX RI FL
        """

rule mapping_ann:
    input:
        fq = "03_pychopper/{sample}_full_length.fq.gz",
        tfa = "10_gffread/gff_cmp.annotated_transcriptome.fa"
    output:
        osam = "12_ann_mapping/{sample}.sam",
        obam = "12_ann_mapping/{sample}.bam",
        obai = "12_ann_mapping/{sample}.bam.bai",
        tse = "13_tse/{sample}.tsv"
    conda:
        "Envs/minimap2.yaml"
    threads: 12
    shell:
        """
        minimap2 -t {threads}  -ax map-ont -p 0 -N 10 {input.tfa} {input.fq} > {output.osam}
        samtools view -@ {threads} -bh {output.osam} | samtools sort -@ {threads} -O BAM - > {output.obam}
        samtools index {output.obam}
        /home/dguan/.local/bin/NanoCount -i {output.obam} -o {output.tse} -3 50 --extra_tx_info
        """

rule run_cppred:
    input:
        tfa = "10_gffread/gff_cmp.annotated_transcriptome.fa",
    output:
        cppred = "14_cppred/gff_cmp.annotated_transcriptome.cppred.result"
    params:
        path="99_scripts/CPPred/bin/"
    conda:
        "Envs/cppred.yaml"
    threads: 1
    shell:
        """
        cp {input.tfa} {params.path}/tmp.fa
        cd {params.path}
        python2 CPPred.py -i tmp.fa -hex ../Hexamer/Human_Hexamer.tsv -r ../Human_Model/Human.range -mol ../Human_Model/Human.model -spe Human -o tmp.cppred.result
        cp tmp.cppred.result {output.cppred}
        """
