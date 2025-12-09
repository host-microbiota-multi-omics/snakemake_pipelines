import os

configfile: "5_metagenomics3.yaml"

# Config variables
WORKDIR = config["workdir"]
BINDIR = config["bins"]

# List bins to process
MAGS, = glob_wildcards(f"{BINDIR}/{{mag}}.fa")

rule all:
    input:
        f"{WORKDIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
        expand(f"{WORKDIR}/metagenomics/kegg/{{mag}}.tsv", mag=MAGS)

rule gtdbtk:
    output:
        f"{WORKDIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        bindir=BINDIR,
        outdir=f"{WORKDIR}/metagenomics/gtdbtk",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 120 * 2 ** (attempt - 1)
    conda:
        "/projects/alberdilab/data/environments/drakkar/afdf1f313a19f553d72dc129bb35c2f8_"
    shell:
        """
        export GTDBTK_DATA_PATH=/datasets/globe_databases/gtdbtk_db/20241001
        gtdbtk classify_wf \
            --genome_dir {params.bindir} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --extension fa \
            --skip_ani_screen
        """

rule prodigal:
    input:
        f"{BINDIR}/{{mag}}.fa"
    output:
        gff=f"{WORKDIR}/metagenomics/prodigal/{{mag}}.gff",
        nt=f"{WORKDIR}/metagenomics/prodigal/{{mag}}.fna",
        aa=f"{WORKDIR}/metagenomics/prodigal/{{mag}}.faa"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb * 10) * 2 ** (attempt - 1))
    shell:
        """module 
        mkdir -p $(dirname {output.gff})
        module load pprodigal/1.0.1
        prodigal -i {input} -o {output.gff} -d {output.nt} -a {output.aa} -p single
        """

rule kegg:
    input:
        f"{WORKDIR}/metagenomics/prodigal/{{mag}}.faa"
    output:
        txt=f"{WORKDIR}/metagenomics/kegg/{{mag}}.txt",
        tsv=f"{WORKDIR}/metagenomics/kegg/{{mag}}.tsv"
    params:
        db="/projects/alberdilab/data/databases/drakkar/kofams"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 1024 * 4) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(20, int(input.size_mb * 60) * 2 ** (attempt - 1))
    threads: 1
    shell:
        """
        module load hmmer/3.3.2
        hmmscan -o {output.txt} --tblout {output.tsv} -E 1e-10 --noali {params.db} {input}
        """