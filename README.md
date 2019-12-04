# Tiger Salamander Project
SNP calling, haplotyping, and subsampling pipeline using [Snakemake](http://snakemake.readthedocs.io/en/stable/index.html) for amplicon sequence data.
Written for the Tiger Salamander project in the [Weisrock Lab](http://sweb.uky.edu/~dweis2/The_Weisrock_Lab/Front_Page.html) at the University of Kentucky.

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.11.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/kelly-sovacool/tiger_salamander_project/blob/master/LICENSE.txt)

## Setup

You can download this repository with:
```
$ git clone https://github.com/kelly-sovacool/tiger_salamander_project.git
```

I recommend using the [Conda](https://conda.io/docs/) package manager. If you don't already have Conda installed, the fastest way to get up and running is to use the [Miniconda](https://conda.io/miniconda.html) Python 3 distribution, which includes Conda.

After installing Conda, change into the directory containing the repository. Then, create a Conda environment with:
```
$ conda env create --name tiger_salamander_project --file config/environment.yml
```

Conda will create an environment with all the dependencies specified in the environment file. To activate the environment, run:
```
$ source activate tiger_salamander_project
```

Alternatively, if you would prefer to use a different package management tool (e.g. `pip3`) and/or install packages system-wide, you can manually install the dependencies listed in `config/environment.yml`.

## Usage

There are two main pipelines here: `haplotype_pipeline` and `snp_pipeline`. The general workflow is to call variants and haplotypes with the `haplotype_pipeline`, manually curate the alignments (optional), then call SNPs and create subsamples with the `snp_pipeline` before feeding the results into a program such as Structure.
The project directory structure is as follows:
```
.
├── haplotype_pipeline
│   ├── adapters
│   ├── reports
│   └── rules
├── legacy_scripts
├── reference
└── snp_pipeline
    ├── haplotypes_curated
    ├── reports
    └── scripts
```

### Haplotype Pipeline

If you use the default config.yaml file, place the Illumina fastq files in `haplotype_pipeline/data/illumina_fastq`, with a separate subdirectory for each sequencing run. Place the 454 haplotype fasta files in `haplotype_pipeline/data/454_haplotypes`.
Finally, change into the `haplotype_pipeline` directory and run it using as many cores as are available with:
```
$ snakemake -j
```

The haplotype pipeline will output haplotypes as single-locus fasta files in `haplotype_pipeline/haplotypes`.
If desired, these can be manually curated with a tool such as Geneious.

### SNP Pipeline

Place curated haplotypes in `snp_pipeline/haplotypes_curated` (alternatively, copy the files from `haplotype_pipeline/haplotypes`). Be sure to edit `config.yaml` so that it contains your email address and the absolute path to Structure on your DLX account. Run it with the same `snakemake` command as above from the `snp_pipeline` directory.

The `snp_pipeline` outputs filtered SNP sites as single-locus fasta files in `snp_pipeline/snp_sites_filtered`, subsamples of the SNP data in Structure format and scripts for running Structure on the DLX in `snp_pipeline/snp_subsamples`. The subsamples and scripts can be copied to your DLX account for running Structure with:
```
$ scp -r snp_pipeline/snp_subsamples username@server:~/
```

### Report

In the top-level directory is a Snakefile with a single rule to generate a report as an html file. Run `snakemake` from the top-level directory to generate it. This report summarizes the results of both pipelines and contains links to histograms for visualizing individuals in loci and SNP sites per locus. Examples:

![alt text](https://github.com/kelly-sovacool/tiger_salamander_project/blob/master/haplotype_pipeline/reports/loci_histogram.png)
![alt text](https://github.com/kelly-sovacool/tiger_salamander_project/blob/master/snp_pipeline/reports/snp_histogram.png)

### Notes:
 * This pipeline was written for a SNP-only analysis. The `filter_variants` rule in `rules/haplotype_illumina_data.smk` will filter out indels. If you wish to keep indels, you can remove the `--remove-indels` flag in that rule.
 * Various scripts and small programs from before the time of Snakemake are now in `legacy_scripts/` for posterity.
