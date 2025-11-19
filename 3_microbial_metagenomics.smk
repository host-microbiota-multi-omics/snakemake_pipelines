import os

configfile: "3_microbial_metagenomics.yaml"

# Config variables
WORKDIR = config["workdir"]
READS = config["reads"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{READS}/{{sample}}_1.fq.gz")

rule all:
    input:
        expand(f"{WORKDIR}/metagenomics/metabat2/{{sample}}.tsv", sample=SAMPLES),
        expand(f"{WORKDIR}/metagenomics/maxbin2/{{sample}}.summary", sample=SAMPLES)

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
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.sample} using metabat2..."
    shell:
        """
        module load metabat2/2.17
        metabat2 -i {input.assembly} -a {input.depth} -o {output} -m 1500 --saveCls --noBinOut
        """

rule maxbin2:
    input:
        assembly=f"{WORKDIR}/metagenomics/megahit/{{sample}}.fna",
        depth=f"{WORKDIR}/metagenomics/bowtie2/{{sample}}_maxbin.depth"
    output:
         f"{WORKDIR}/metagenomics/maxbin2/{{sample}}.summary"
    params:
        basename=f"{WORKDIR}/metagenomics/maxbin2/{{sample}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 3) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.sample} using maxbin2..."
    shell:
        """
        module load maxbin2/2.2.7 hmmer/3.3.2
        rm -rf {params.basename}*
        run_MaxBin.pl -contig {input.assembly} -abund {input.depth} -max_iteration 10 -out {params.basename} -min_contig_length 1500
        """

rule semibin2:
    input:
        assembly=f"{WORKDIR}/metagenomics/megahit/{{sample}}.fna",
        bam=f"{WORKDIR}/metagenomics/bowtie2/{{sample}}.bam"
    output:
        f"{WORKDIR}/metagenomics/semibin2/{{sample}}/contig_bins.tsv"
    params:
        outdir=f"{WORKDIR}/metagenomics/semibin2/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: min(1000*1024,max(8*1024, int(input.size_mb * 30) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb / 2) * 2 ** (attempt - 1)))
    message: "Binning contigs from assembly {wildcards.sample} using semibin2..."
    shell:
        """
        module load {params.semibin2_module} {params.bedtools_module} {params.hmmer_module}
        SemiBin2 single_easy_bin -i {input.assembly} -b {input.bam} -o {params.outdir} -m 1500 -t {threads} --compression none
        """

rule semibin2_table:
    input:
         f"{WORKDIR}/metagenomics/semibin2/{{sample}}/contig_bins.tsv"
    output:
        f"{WORKDIR}/metagenomics/semibin2/{{sample}}.tsv"
    params:
        fastadir=f"{WORKDIR}/metagenomics/semibin2/{{sample}}/output_bins"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(5, int(input.size_mb / 5) * 2 ** (attempt - 1))
    shell:
        """
        tail -n +2 {input} > {output}
        """


# Script to calculate resources based on the number of bins
_row_count_cache = {}
def row_count(path):
    """Return number of data rows (excluding header) in a TSV, caching the result."""
    if path not in _row_count_cache:
        with open(path) as f:
            _row_count_cache[path] = max(0, sum(1 for _ in f))
    return _row_count_cache[path]

checkpoint binette:
    input:
        metabat2=f"{WORKDIR}/metagenomics/metabat2/{{sample}}.tsv",,
        maxbin2=f"{WORKDIR}/metagenomics/maxbin2/{{sample}}.tsv",
        semibin2=f"{WORKDIR}/metagenomics/semibin2/{{sample}}.tsv",
        fasta=f"{WORKDIR}/metagenomics/megahit/{{sample}}.fna"
    output:
        f"{WORKDIR}/metagenomics/binette/final_bins_quality_reports.tsv"
    params:
        checkm_db = {CHECKM2_DB},
        diamond_module = {DIAMOND_MODULE},
        checkm2_module = {CHECKM2_MODULE},
        binette_module = {BINETTE_MODULE},
        outdir=f"{OUTPUT_DIR}/cataloging/binette/{{sample}}"
    threads: 8
    conda:
        f"{PACKAGE_DIR}/workflow/envs/cataloging.yaml"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(1000*1024,max(32*1024, (row_count(input.metabat2) + row_count(input.maxbin2) + row_count(input.semibin2)) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb) * 2 ** (attempt - 1)))
    message: "Refining bins from assembly {wildcards.assembly} using binette..."
    shell:
        """
        # Define input files
        METABAT2="{input.metabat2}"
        MAXBIN2="{input.maxbin2}"
        SEMIBIN2="{input.semibin2}"

        # Remove empty input files from the list
        VALID_TSV_FILES=""
        if [ -s "$METABAT2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $METABAT2"
        fi
        if [ -s "$MAXBIN2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $MAXBIN2"
        fi
        if [ -s "$SEMIBIN2" ]; then
            VALID_TSV_FILES="$VALID_TSV_FILES $SEMIBIN2"
        fi

        # Ensure at least one valid TSV file exists
        if [ -z "$VALID_TSV_FILES" ]; then
            echo "Error: No valid TSV input files for binette." >&2
            exit 1
        fi

        # Run binette only with non-empty TSV files
        binette --contig2bin_tables $VALID_TSV_FILES \
                --contigs {input.fasta} \
                --outdir {params.outdir} \
                --checkm2_db {params.checkm_db} \
                --threads {threads}
        """

# Regenerate the bin_id wildcard based on the checkpoint results
def get_bin_fna_sep(wildcards):
    checkpoint_output = checkpoints.binette.get(**wildcards).output[0]
    cluster_ids = get_bin_ids_from_tsv(checkpoint_output)
    return f"{OUTPUT_DIR}/cataloging/binette/{{sample}}/final_bins/bin_{wildcards.bin_id}.fa"

rule rename_bins:
    input:
        lambda wildcards: get_bin_fna_sep(wildcards)
    output:
        f"{OUTPUT_DIR}/cataloging/binette/{{sample}}/final_bins/bin_{{bin_id}}.fa"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(2, int(input.size_mb) * 2 ** (attempt - 1))
    message: "Renaming bin {wildcards.bin_id} from assembly {wildcards.assembly}..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/rename_bins.py {wildcards.assembly} {input} {output}
        """

rule move_metadata:
    input:
        f"{OUTPUT_DIR}/cataloging/binette/{{assembly}}/final_bins_quality_reports.tsv"
    output:
        f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv"
    threads: 1
    resources:
        mem_mb=8*1024,
        runtime=10
    message: "Exporting bin metadata from assembly {wildcards.assembly}..."
    shell:
        """
        cp {input} {output}
        """

rule all_bins:
    input:
        expand(f"{OUTPUT_DIR}/cataloging/final/{{assembly}}.tsv", assembly=assemblies)
    output:
        paths=f"{OUTPUT_DIR}/cataloging/final/all_bin_paths.txt",
        metadata=f"{OUTPUT_DIR}/cataloging/final/all_bin_metadata.csv"
    params:
        package_dir={PACKAGE_DIR}
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(2, int(input.size_mb) * 2 ** (attempt - 1))
    message: "Generating bin path file..."
    shell:
        """
        python {params.package_dir}/workflow/scripts/all_bin_paths.py {input} -o {output.paths}
        python {params.package_dir}/workflow/scripts/all_bin_metadata.py {input} -o {output.metadata}
        """
