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
