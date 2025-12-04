import os

configfile: "2_genomics.yaml"

# Config variables
WORKDIR = config["workdir"]
BAMS = config["bams"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{BAMS}/{{sample}}.bam")

rule all:
    input:
        f"{WORKDIR}/genomics/pcoa.png",
        f"{WORKDIR}/genomics/heatmap.png"

rule bam_to_fastq:
    input:
        f"{BAMS}/{{sample}}.bam"
    output:
        r1=f"{WORKDIR}/genomics/reads/{{sample}}_1.fq",
        r2=f"{WORKDIR}/genomics/reads/{{sample}}_2.fq"
    threads: 4
    params:
        collate="",
        fastq="-n",
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    shell:
        """
        module load samtools/1.21
        samtools collate -u -O {input} | samtools fastq -1 {output.r1} -2 {output.r2}
        """

rule simplify_reads:
    input:
        r1=f"{WORKDIR}/genomics/reads/{{sample}}_1.fq",
        r2=f"{WORKDIR}/genomics/reads/{{sample}}_2.fq"
    output:
        f"{WORKDIR}/genomics/reads/{{sample}}/{{sample}}.fq"
    threads: 1
    localrule: True
    shell:
        """
        rm {input.r2}
        mv {input.r1} {output}
        """

rule skmer_reference:
    input:
        f"{WORKDIR}/genomics/reads/{{sample}}/{{sample}}.fq"
    output:
        f"{WORKDIR}/genomics/skmer/{{sample}}/{{sample}}.dat"
    threads: 1
    params:
        inputdir=f"{WORKDIR}/genomics/reads/{{sample}}",
        outputbase="skmer"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 60) * 2 ** (attempt - 1))
    shell:
        """
        module load skmer/3.3.0
        skmer reference {params.inputdir} -p {threads} -l {params.outputbase}
        """

rule skmer_distance:
    input:
        expand(f"{WORKDIR}/genomics/skmer/{{sample}}/{{sample}}.dat", sample=SAMPLES)
    output:
        f"{WORKDIR}/genomics/distance_matrix.txt"
    threads: 1
    params:
        inputdir=f"{WORKDIR}/genomics/skmer",
        outputbase="distance_matrix"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 50) * 2 ** (attempt - 1))
    shell:
        """
        module load skmer/3.3.0
        skmer distance {params.inputdir} -p {threads} -o {params.outputbase}
        """

rule skmer_pcoa:
    input:
        f"{WORKDIR}/genomics/distance_matrix.txt"
    output:
        f"{WORKDIR}/genomics/pcoa.png"
    log:
        "logs/skmer/pcoa.log"
    threads: 1
    conda:
        "envs/plot.yml"
    params:
        outputbase="distance_matrix"
    resources:
        mem_mb=lambda wildcards, input, attempt: 4*1024,
        runtime=lambda wildcards, input, attempt: 5
    shell:
        """
        python workflow/scripts/pcoa_from_distance.py --dist {input} --out-prefix pcoa
        """

rule skmer_heatmap:
    input:
        f"{WORKDIR}/genomics/distance_matrix.txt"
    output:
        f"{WORKDIR}/genomics/heatmap.png"
    threads: 1
    conda:
        "envs/plot.yml"
    params:
        outputbase="distance_matrix"
    resources:
        mem_mb=lambda wildcards, input, attempt: 4*1024,
        runtime=lambda wildcards, input, attempt: 5
    shell:
        """
        python workflow/scripts/heatmap_from_distance.py --dist {input} --out {output}
        """
