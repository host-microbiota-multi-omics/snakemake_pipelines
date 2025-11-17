configfile: "1_preprocessing.yaml"

# Config variables
WORKDIR = config["workdir"]
READS = config["reads"]
REFERENCE = config["reference"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{READS}/{{sample}}_1.fq.gz")

rule all:
    input:
        expand(f"{WORKDIR}/preprocessing/{{sample}}_1.fq.gz", sample=SAMPLES),
        expand(f"{WORKDIR}/preprocessing/{{sample}}_2.fq.gz", sample=SAMPLES),
        expand(f"{WORKDIR}/preprocessing/{{sample}}.bam", sample=SAMPLES)

rule fastp:
    input:
        r1=f"{READS}/{{sample}}_1.fq.gz",
        r2=f"{READS}/{{sample}}_2.fq.gz"
    output:
        r1=f"{WORKDIR}/preprocessing/fastp/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/preprocessing/fastp/{{sample}}_2.fq.gz",
        html=f"{WORKDIR}/preprocessing/fastp/{{sample}}.html",
        json=f"{WORKDIR}/preprocessing/fastp/{{sample}}.json"
    params:
        q=config["qualified_quality_phred"],
        l=config["length_required"]
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 3) * 2 ** (attempt - 1))
    message: "Quality-filtering sample {wildcards.sample}..."
    shell:
        """
        module load fastp/0.23.4
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --trim_poly_g \
            --trim_poly_x \
            --low_complexity_filter \
            --n_base_limit 5 \
            --qualified_quality_phred {params.q} \
            --length_required {params.l} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        """

# Build a Bowtie2 index from the given reference genome FASTA file.  
# Produces the index files required for read mapping and saves a copy of the reference sequence.

rule reference_index:
    input:
        REFERENCE
    output:
        index=f"{WORKDIR}/reference/{{reference}}.rev.1.bt2"
    params:
        basename=f"{WORKDIR}/reference/{{reference}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Indexing reference genome {wildcards.reference}..."
    shell:
        """
        module load bowtie2/2.4.2
        bowtie2-build {input} {params.basename}
        cat {input} > {params.basename}.fna
        """

# Align quality-filtered paired-end reads to the corresponding reference genome using Bowtie2.  
# The alignments are converted to BAM format and sorted with samtools for downstream analyses.

rule reference_map:
    input:
        index=lambda wildcards: expand(
            f"{WORKDIR}/reference/{{reference}}.rev.1.bt2",
            reference=[SAMPLE_TO_REFERENCE[wildcards.sample]]
        ),
        r1=f"{WORKDIR}/preprocessing/fastp/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/preprocessing/fastp/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/preprocessing/bowtie2/{{sample}}.bam"
    params:
        basename=lambda wildcards: f"{WORKDIR}/reference/{SAMPLE_TO_REFERENCE[wildcards.sample]}"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against reference genome..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """

# Generate alignment quality metrics for each BAM file using samtools.  
# Produces an index (.bai), a flagstat summary, idxstats, and detailed stats file for MultiQC reporting.

rule samtools_stats:
    input:
        rules.reference_map.output
    output:
        bai      = f"{WORKDIR}/preprocessing/samtools/{{sample}}.bam.bai",
        flagstat = f"{WORKDIR}/preprocessing/samtools/{{sample}}.flagstat.txt",
        idxstats = f"{WORKDIR}/preprocessing/samtools/{{sample}}.idxstats.txt",
        stats    = f"{WORKDIR}/preprocessing/samtools/{{sample}}.stats.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 2) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: 10 * 2 ** (attempt - 1)
    message: "Generating mapping stats for {wildcards.sample}..."
    shell:
        """
        module load samtools/1.21
        samtools index {input} {output.bai}
        samtools flagstat {input} > {output.flagstat}
        samtools idxstats {input} > {output.idxstats}
        samtools stats {input} > {output.stats}
        """

# Split mapped BAM files into metagenomic (unmapped) and host (mapped) read sets.  
# Outputs paired FASTQ files for unmapped reads, a host-only BAM, and text files counting reads/bases in each category.

rule split_reads:
    input:
        f"{WORKDIR}/preprocessing/bowtie2/{{sample}}.bam"
    output:
        r1=f"{WORKDIR}/preprocessing/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/preprocessing/{{sample}}_2.fq.gz",
        bam=f"{WORKDIR}/preprocessing/{{sample}}.bam"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Extracting metagenomic reads of {wildcards.sample}..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        samtools view -b -f12 -@ {threads} {input} | samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -
        samtools view -b -f12 -@ {threads} {input} | samtools view -c - > {output.metareads}
        samtools view -f12 -@ {threads} {input} | awk '{{sum += length($10)}} END {{print sum}}' > {output.metabases}
        samtools view -b -F12 -@ {threads} {input} | samtools sort -@ {threads} -o {output.bam} -
        samtools view -b -F12 -@ {threads} {input} | samtools view -c - > {output.hostreads}
        samtools view -F12 -@ {threads} {input} | awk '{{sum += length($10)}} END {{print sum}}' > {output.hostbases}
        """

rule singlem:
    input: 
        r1=f"{WORKDIR}/preprocessing/singlem/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/preprocessing/singlem/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/preprocessing/singlem/{{sample}}.profile"
    params:
        siglemdir = f"{WORKDIR}/preprocessing/singlem/"
    threads: 1
    shell:
        """
        module load singlem/0.19.0
        export SINGLEM_METAPACKAGE_PATH=/maps/datasets/globe_databases/singlem/5.4.0/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
        mkdir -p {params.siglemdir}
        rm -rf {params.workdir}
        singlem pipe \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {output}
        """

rule spf:
    input: 
        r1=f"{WORKDIR}/preprocessing/singlem/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/preprocessing/singlem/{{sample}}_2.fq.gz",
        profile=f"{WORKDIR}/preprocessing/singlem/{{sample}}.profile"
    output:
        f"{WORKDIR}/preprocessing/singlem/{{sample}}.fraction"
    params:
        workdir = lambda wc: f"{WORKDIR}/preprocessing/singlem/{wc.sample}"
    threads: 1
    shell:
        """
        module load singlem/0.19.0
        export SINGLEM_METAPACKAGE_PATH=/maps/datasets/globe_databases/singlem/5.4.0/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
        singlem microbial_fraction \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {input.profile} > {output}
        """
