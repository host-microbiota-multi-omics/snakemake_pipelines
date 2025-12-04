import os

configfile: "4_annotations.yaml"

# Config variables
WORKDIR = config["workdir"]
BINS = config["bins"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{READS}/{{sample}}_1.fq.gz")

rule all:
    input:
        f"{OUTPUT_DIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv"

rule gtdbtk:
    output:
        f"{OUTPUT_DIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        bindir=BINS,
        outdir=f"{WORKDIR}/metagenomics/gtdbtk",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 120 * 2 ** (attempt - 1)
    shell:
        """
        export GTDBTK_DATA_PATH=/datasets/globe_databases/gtdbtk_db/20241001
        gtdbtk classify_wf \
            --genome_dir {params.bindir} \
            --out_dir {params.outdir} \
            --cpus {threads}
        """