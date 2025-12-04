import os

configfile: "3_metagenomics.yaml"

# Config variables
WORKDIR = config["workdir"]
READS = config["reads"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{READS}/{{sample}}_1.fq.gz")

rule all:
    input:
        expand(f"{WORKDIR}/metagenomics/drep/{{sample}}/data_tables/genomeInformation.csv", sample=SAMPLES)

rule assembly:
    input:
        r1=f"{READS}/{{sample}}_1.fq.gz",
        r2=f"{READS}/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/metagenomics/megahit/{{sample}}.fna"
    params:
        outputdir=f"{WORKDIR}/metagenomics/megahit/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: min(1020*1024,max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Assembling {wildcards.sample}..."
    shell:
        """
        module load megahit/1.2.9
        rm -rf {params.outputdir}

        megahit \
            -t {threads} \
            --verbose \
            --min-contig-len 1500 \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outputdir}
        mv {params.outputdir}/final.contigs.fa {output}
        """

rule assembly_index:
    input:
        f"{WORKDIR}/metagenomics/megahit/{{sample}}.fna"
    output:
        index=f"{WORKDIR}/metagenomics/megahit/{{sample}}.rev.2.bt2"
    params:
        basename=f"{WORKDIR}/metagenomics/megahit/{{sample}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Indexing assembly {wildcards.sample}..."
    shell:
        """
        module load bowtie2/2.4.2
        bowtie2-build {input} {params.basename}
        """

rule assembly_map:
    input:
        index=f"{WORKDIR}/metagenomics/megahit/{{sample}}.rev.2.bt2",
        r1=f"{READS}/{{sample}}_1.fq.gz",
        r2=f"{READS}/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/metagenomics/bowtie2/{{sample}}.bam"
    params:
        basename=f"{WORKDIR}/metagenomics/megahit/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 3) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} reads to assembly..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """

rule assembly_map_depth:
    input:
        f"{WORKDIR}/metagenomics/bowtie2/{{sample}}.bam"
    output:
        metabat2=f"{WORKDIR}/metagenomics/bowtie2/{{sample}}_metabat.depth",
        maxbin2=f"{WORKDIR}/metagenomics/bowtie2/{{sample}}_maxbin.depth"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Calculating mapping states of assembly {wildcards.sample}..."
    shell:
        """
        module load metabat2/2.17
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """

rule metabat2:
    input:
        assembly=f"{WORKDIR}/metagenomics/megahit/{{sample}}.fna",
        depth=f"{WORKDIR}/metagenomics/bowtie2/{{sample}}_metabat.depth"
    output:
        f"{WORKDIR}/metagenomics/metabat2/{{sample}}.tsv"
    params:
        basename=f"{WORKDIR}/metagenomics/metabat2/{{sample}}/{{sample}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.sample} using metabat2..."
    shell:
        """
        module load metabat2/2.17
        metabat2 -i {input.assembly} -a {input.depth} -o {params.basename} -m 1500 --saveCls

        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "*$(basename {params.basename})_*.fa" | sort > {output}
        """

rule checkm:
    input:
        f"{WORKDIR}/metagenomics/metabat2/{{sample}}.tsv"
    output:
        f"{WORKDIR}/metagenomics/checkm2/{{sample}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/metagenomics/metabat2/{wildcards.sample}",
        outdir=f"{WORKDIR}/metagenomics/checkm2/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fa -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "genome","completeness","contamination"; next}} {{print $1".fna",$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule drep:
    input:
        genomes=f"{WORKDIR}/metagenomics/metabat2/{{sample}}.tsv",
        genomeinfo=f"{WORKDIR}/metagenomics/checkm2/{{sample}}.tsv"
    output:
        f"{WORKDIR}/metagenomics/drep/{{sample}}/data_tables/genomeInformation.csv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/metagenomics/metabat2/{wildcards.sample}",
        outdir=f"{WORKDIR}/metagenomics/drep/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3
        rm -rf {params.outdir}
        dRep dereplicate {params.outdir} -g {input.genomes} -p {threads} -pa 0.95 --genomeInfo {input.genomeinfo}
        """