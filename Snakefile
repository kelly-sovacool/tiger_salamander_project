""" SNP calling, haplotyping, and subsampling pipeline for the Tiger Salamander Project """
import collections
import HTSeq
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
#from xhtml2pdf import pisa

subworkflow haps:
    workdir: "haplotype_pipeline"

subworkflow snps:
    workdir: "snp_pipeline"

wildcards = glob_wildcards(haps("intermediates_illumina/alignments/{sample}.sorted.bam"))
samples = set(wildcards.sample)

rule targets:
    input:
        read_counts_filtered=snps("results/read_counts_filtered.csv"),
        haps_html=haps("results/loci_histogram.html"),
        haps_mtx=haps("results/matrix.tsv"),
        snps_html=snps("results/snp_histogram.html"),
        snps_log=snps("logs/histogram/snps_per_locus.log")