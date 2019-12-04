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

rule html_to_pdf:  # TODO: convert plotly html files to pdf. or use seaborn intead of plotly to create SVGs
    input:
        haps_hist=haps("results/loci_histogram.html"),
        snps_hist=snps("results/snp_histogram.html")
    output:
        haps_pdf="results/loci_histogram.pdf",
        snps_pdf="results/snp_histogram.pdf"
    run:  # TODO: xhtml2pdf requires python 2.7
        for input_fn, output_fn in ((input.haps_hist, output.haps_pdf), (input.snps_hist, output.snps_hist)):
            with open(input_fn, 'r') as input_file:
                source_html = input_file.read()
            with open(output_fn, "w+b") as output_file:
                pisa_status = pisa.CreatePDF(source_html, dest=output_file)
            # return True on success and False on errors
            if not pisa_status.err:
                raise ValueError("pisa status error {}".format(pisa_status.err))

rule report:
    input:
        haps_pdf="results/loci_histogram.pdf",
        haps_mtx=haps("results/matrix.tsv"),
        snps_pdf="results/snp_histogram.pdf",
        snps_log=snps("logs/histogram/snps_per_locus.log")
    output:
        report="results/report.html"
    run:
        from snakemake.utils import report
        with open(input.haps_mtx) as file:
            matrix = file.readlines()
        num_individuals = len(matrix) - 1
        num_loci = len(matrix[0]) - 2
        with open(input.snps_log, 'r') as file:
            avg_snps = round(float(file.read().strip()),2)
        report("""
        Pipeline for the Tiger Salamander Project
        =========================================
        Includes SNP calling, haplotyping, and subsampling.

        {num_individuals} individuals were haplotyped for {num_loci} (see Graph G1_ and Table T1_).

        Average of {avg_snps} SNPs per locus (see Graph G2_).

        {loci_hist}

        {snp_hist}
        """, output.report, T1=input.haps_mtx, G1=input.haps_pdf, G2=input.snps_pdf)
