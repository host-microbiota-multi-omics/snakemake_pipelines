import os

configfile: "4_metagenomics2.yaml"

# Config variables
WORKDIR = config["workdir"]
BINS = config["bins"]
READS = config["reads"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{READS}/{{sample}}_1.fq.gz")

rule all:
    input:
        #f"{WORKDIR}/metagenomics/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
        f"{WORKDIR}/metagenomics/coverm/coverm_read_counts.tsv",
        f"{WORKDIR}/metagenomics/coverm/coverm_covered_bases.tsv"

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
        export GTDBTK_DATA_PATH=/datasets/globe_databases/gtdbtk_db/20241001
        gtdbtk classify_wf \
            --genome_dir {params.bindir} \
            --out_dir {params.outdir} \
            --cpus {threads}
        """

rule concatenate_catalogue:
    output:
        f"{WORKDIR}/metagenomics/catalogue/catalogue.fna"
    params:
        bindir=BINS
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 128*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 120 * 2 ** (attempt - 1)
    shell:
        """
        cat {params.bindir}/*.fna > {output}
        """

rule index_catalogue:
    input:
        f"{WORKDIR}/metagenomics/catalogue/catalogue.fna"
    output:
        f"{WORKDIR}/metagenomics/catalogue/catalogue.fna.bt2"
    params:
        basename = f"{WORKDIR}/metagenomics/catalogue/catalogue"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 64*1024 * 2 ** (attempt - 1),
        runtime=lambda wildcards, attempt: 60 * 2 ** (attempt - 1)
    shell:
        """
        module load bowtie2/2.4.2
        bowtie2-build {input} {params.basename}
        """

rule map_reads_to_catalogue:
    input:
        r1=f"{WORKDIR}/preprocessing/final/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/preprocessing/final/{{sample}}_2.fq.gz",
        index=f"{WORKDIR}/metagenomics/catalogue/catalogue.fna.bt2"
    output:
        bam=f"{WORKDIR}/metagenomics/mapping/{{sample}}.bam"
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(16*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(30, int(input.size_mb / 100) * 2 ** (attempt - 1))
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        bowtie2 -x {input.index} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -b -@ {threads} -o {output.bam} -
        """

rule quantify_reads_catalogue:
    input:
        expand(f"{WORKDIR}/metagenomics/mapping/{{sample}}.bam", sample=samples)
    output:
        f"{WORKDIR}/metagenomics/coverm/coverm.tsv"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb / 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 / 5) * 2 ** (attempt - 1))
    message:
        "Generating mapping statistics with..."
    shell:
        """
        mkdir -p $(dirname {output})
        module load coverm/0.7.0
        coverm genome \
            -b {input} \
            -s ^ \
            -m count covered_bases \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}
        """

rule split_coverm:
    input:
        coverm=f"{WORKDIR}/metagenomics/coverm/coverm.tsv"
    output:
        read_counts=f"{WORKDIR}/metagenomics/coverm/coverm_read_counts.tsv",
        covered_bases=f"{WORKDIR}/metagenomics/coverm/coverm_covered_bases.tsv"
    localrule: True
    run:
        from pathlib import Path
        import pandas as pd

        read_counts_path = Path(output.read_counts)
        covered_bases_path = Path(output.covered_bases)
        read_counts_path.parent.mkdir(parents=True, exist_ok=True)
        covered_bases_path.parent.mkdir(parents=True, exist_ok=True)

        df = pd.read_csv(input.coverm, sep="\t")
        genome_col = df.iloc[:, 0]

        read_counts = {"Genome": genome_col}
        covered_bases = {"Genome": genome_col}

        for col in df.columns[1:]:
            if "Read Count" in col:
                read_counts[col.replace(" Read Count", "")] = df[col]
            elif "Covered Bases" in col:
                covered_bases[col.replace(" Covered Bases", "")] = df[col]

        pd.DataFrame(read_counts).to_csv(read_counts_path, sep="\t", index=False)
        pd.DataFrame(covered_bases).to_csv(covered_bases_path, sep="\t", index=False)
  
