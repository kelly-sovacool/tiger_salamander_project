""" SNP calling, haplotyping, and subsampling pipeline for the Tiger Salamander Project """
import os

subworkflow haps:
    workdir: "haplotype_pipeline"

subworkflow snps:
    workdir: "snp_pipeline"

rule report:
    input:
        haps_hist=haps("reports/loci_histogram.html"),
        haps_mtx=haps("reports/matrix.tsv")
        snps_hist=snps("reports/snp_histogram.html"),
        snps_log=snps("logs/histogram/snps_per_locus.log")
    output:
        report="reports/report.html"
    run:
        from snakemake.utils import report
        with open(input.haps_mtx) as file:
            matrix = file.readlines()
        num_individuals = len(matrix) - 1
        num_loci = len(matrix[0]) - 2
        with open(input.snps_log, 'r') as file:
            avg_snps = file.read().strip()
        report("""
        SNP calling, haplotyping, and subsampling pipeline for the Tiger Salamander Project
        ===================================
        {num_individuals} individuals were haplotyped for {num_loci} loci (see Table T1_ and Graph G1_).
        Average of {avg_snps} SNPs per locus (see Graph G2_).
        """, output.report, T1=input.haps_mtx, G1=input.haps_hist, G2=input.snps_hist)
