# Snakemake pipelines

This repository contains a collection of snakemake pipelines optimised for Mjolnir.

## Get code

```sh
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/config.yaml
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.smk
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.yaml
```

```sh
module load snakemake/9.9.0
snakemake -s {WORKDIR} --cores 8 --profile config.yaml
```
