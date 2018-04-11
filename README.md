# Tiger Salamander Project
SNP calling, haplotype, and subsampling pipeline using [Snakemake](http://snakemake.readthedocs.io/en/stable/index.html) for amplicon sequence data.
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

If you use the default config.yaml file, place the Illumina fastq files in `haplotype_pipeline/data/illumina_fastq`, with a separate subdirectory for each sequencing run. Place the 454 haplotype fasta files in `haplotype_pipeline/data/454_haplotypes`.
Finally, change into the `haplotype_pipeline` directory and run it using as many cores as are available with:
```
$ snakemake -j
```

The haplotype pipeline will output haplotypes as single-locus fasta files in `haplotype_pipeline/haplotypes`.
If desired, these can be manually curated with a tool such as Geneious.

Place curated haplotypes in `snp_pipeline/haplotypes_curated` (or copy the files from `haplotype_pipeline/haplotypes`) and change into the `snp_pipeline` directory. Be sure to edit `snp_pipeline/config.yaml` so that it contains your email address and path to Structure on your DLX account. Run it with the same command as above.

The `snp_pipeline` outputs filtered SNP sites as single-locus fasta files in `snp_pipeline/snp_sites_filtered`, subsamples of the SNP data in Structure format and scripts for running Structure on the DLX in `snp_pipeline/snp_subsamples`. The subsamples and scripts can be copied to your DLX account for running Structure with:
```
$ scp -r snp_pipeline/snp_subsamples username@server:~/
```

In the top-level directory is a Snakefile with a single rule to generate a report. [This report](report.html) summarizes the results of both pipelines and contains links to histograms for visualizing [individuals in loci](haplotype_pipeline/reports/loci_histogram.html) and [SNP sites per locus](snp_pipeline/reports/snp_histogram.html).


### Notes:
 * This pipeline was written for a SNP-only analysis. The `filter_variants` rule in `rules/haplotype_illumina_data.smk` will filter out indels. If you wish to keep indels, you can remove the `--remove-indels` flag in that rule.
 * Various scripts and small programs from before the time of Snakemake are now in `legacy_scripts/` for posterity.
