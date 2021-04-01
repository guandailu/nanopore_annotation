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
