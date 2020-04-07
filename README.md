# indrops-star
Pipeline to analyze Indrops V3 data using starSolo

*This pipeline is currently under development and should not be considered stable*

## Installation

#### Clone this repository

```git clone https://github.com/BradhamLab/indrops-star/```

#### Install `conda` environment

```conda env create -f environment.yml```

## Run Pipeline

#### Edit Configuration File

The `indrops-star` uses a configuration file `files/config.yaml`. Modify the empty example configuration file `files/empty_config.yaml` in your local clone to suite your needs.

#### Run Pipeline using snakemake

From the head of the repo, execuate the pipeline by issueing the following command:

```
conda activate indrops-star
snakemake
```

For more information for on `Snakemake` see the [documentation](https://snakemake.readthedocs.io/en/stable/) for more information. 
