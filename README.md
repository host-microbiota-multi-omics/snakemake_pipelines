# Snakemake pipelines

This repository contains a collection of snakemake pipelines optimised for Mjolnir.

## Preprocessing

### Get code

```sh
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/config.yaml
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.smk
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.yaml
```

### Prepare the config file

```sh
nano 1_preprocessing.yaml
```

### Execute the pipeline

```sh
screen -S hmmo
module load snakemake/9.9.0
snakemake -s 1_preprocessing.smk --cores 8 --profile .
```
