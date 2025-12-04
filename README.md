# Snakemake pipelines

This repository contains a collection of snakemake pipelines optimised for Mjolnir.

## 0. Prepare the working directory

## Get the Snakemake profile config

```sh
mkdir my_working_dir
cd my_working_dir
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/config.yaml
```

## Create a screen session

Give it `hmmo` (host-microbiota multi-omics) or any other name.

```sh
screen -S hmmo
```

Ctr+a+d to get out from the screen session.

## 1. Preprocessing

### Get the code

```sh
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.smk
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/1_preprocessing.yaml
```

### Prepare the config file

```sh
nano 1_preprocessing.yaml
```

### Execute the pipeline

```sh
screen -r hmmo
module load snakemake/9.9.0
snakemake -s 1_preprocessing.smk --cores 8 --profile .
```

## 2. Host genomics

### Get the code

```sh
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/2_genomics.smk
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/2_genomics.yaml
```

### Prepare the config file

```sh
nano 2_genomics.yaml
```

### Execute the pipeline

```sh
screen -r hmmo
module load snakemake/9.9.0
snakemake -s 2_genomics.smk --cores 8 --profile .
```

## 3. Microbial metagenomics

### Get the code

```sh
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/3_metagenomics.smk
wget https://raw.githubusercontent.com/host-microbiota-multi-omics/snakemake_pipelines/refs/heads/main/3_metagenomics.yaml
```

### Prepare the config file

```sh
nano 3_metagenomics.yaml
```

### Execute the pipeline

```sh
screen -r hmmo
module load snakemake/9.9.0
snakemake -s 3_microbial_metagenomics.smk --cores 8 --profile .
```
