# Tiger Salamander Project
SNP calling & haplotyping pipeline using [Snakemake](http://snakemake.readthedocs.io/en/stable/index.html) for amplicon sequence data.
Written for the Tiger Salamander project in the Weisrock Lab at the University of Kentucky.

## Setup

You can download this repository with:
```
$ git clone https://github.com/kelly-sovacool/tiger_salamander_project.git
```

I recommend using the [Conda](https://conda.io/docs/) package manager. If you don't already have Conda installed, the fastest way to get up and running is to use the [Miniconda](https://conda.io/miniconda.html) Python 3 distribution, which includes Conda.

After installing Conda, change into the directory containing the repository. Then, create a Conda environment with:
```
$ conda env create --name tiger_salamander_project --file environment.yaml
```

Conda will create an environment with all the dependencies specified in the environment.yaml file. To activate the environment, run:
```
$ source activate tiger_salamander_project
```

## Usage

If you use the default config.yaml file, place the Illumina fastq files in `data/illumina_fastq`, with a separate subdirectory for each sequencing run. Place the 454 haplotype fasta files in `data/454_haplotypes`.
Finally, run the pipeline using as many cores as are available with:
```
$ snakemake -j
```

The pipeline will output haplotypes and SNP sites as single-locus fasta files and a table (`reports/matrix.tsv`) showing the presence of samples in each locus.


### Notes:
 * This pipeline was written for a SNP-only analysis. The `filter_variants` rule in `rules/haplotype_illumina_data.smk` will filter out indels. If you wish to keep indels, you can remove the `--remove-indels` flag in that rule.
 * Various scripts and small programs from before the time of Snakemake are now in scripts/ for posterity.
