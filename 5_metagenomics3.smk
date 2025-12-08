import os

configfile: "5_metagenomics3.yaml"

# Config variables
WORKDIR = config["workdir"]
BINS = config["bins"]

rule all:
    input:
        f"{WORKDIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv"

rule gtdbtk:
    output:
        f"{WORKDIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        bindir=BINS,
        outdir=f"{WORKDIR}/metagenomics/gtdbtk",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 120 * 2 ** (attempt - 1)
    shell:
        """
        module load gtdbtk/2.2.6
        export GTDBTK_DATA_PATH=/datasets/globe_databases/gtdbtk_db/20241001
        gtdbtk classify_wf \
            --genome_dir {params.bindir} \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --extension fa \
            --skip_ani_screen
        """
