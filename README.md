# indrops-star
Pipeline to analyze Indrops V3 data using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)

## Installation

#### Clone this repository

```git clone https://github.com/BradhamLab/indrops-star/```

#### Install `conda` environment

```conda env create -f environment.yml```

## Run Pipeline

#### Edit Configuration File

The `indrops-star` uses a configuration file `files/config.yaml`. Modify the empty example configuration file `files/empty_config.yaml` in your local clone to suit your needs. Entries are explained below:

```yaml
genome:
  fasta: "your/genoma/sequences.fasta" # file path to genome fasta file
  gtf: "your/genome/annotations.gtf" # file path to genome annotations as gtf file
STAR:
  index: "star/index/install/dir" # where to build STAR index
project:
  dir: "head/of/your/project/dir/Intensities/Basecalls"  # directory containing fastq files
  id: "your-project-id"  # project identification string
  libraries:  # key value mapping library barcodes to library names
    {"ATTAGAGG": "Library1", # Index was CCTCTAAT, reverse compliment: ATTAGAGG
     "CGGAGAGA": "Library2", # Index was TCTCTCCG <==> CGGAGAGA
     "CTAGTCGA": "Library3", # Index was TCGACTAG <==> CTAGTCGA
     "AGCTAGAA": "Library4"} # Index was TTCTAGCT <==> AGCTAGAA
params:
  weave_fastqs:
    mismatches: 2  # allowable mismatches in library barcodes
  star_index: "--sjdbOverhang 60 --genomeSAIndexNbases 13" # put extra STAR parameters to pass when building index here
  star_solo: "--soloCBmatchWLtype 1MM_multi" # put extra STAR alignment parameters here
```

#### Run Pipeline using snakemake

From the head of the repo, execuate the pipeline by issueing the following command:

```
conda activate indrops-star
snakemake
```

Final STARsolo output will be found in `<your/project/dir>/processed/STAR/`


For more information on `Snakemake` see the [documentation](https://snakemake.readthedocs.io/en/stable/) for more information. 
