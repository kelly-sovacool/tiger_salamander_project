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

rule plot_read_counts:
    input:
        'results/read_counts.csv'
    output:
        heatmap='figures/reads_heatmap.svg',
        barplot_samples='figures/reads_samples.svg',
        barplot_loci='figures/reads_loci.svg'
    run:
        read_counts_df = pd.read_csv(input[0], index_col=0)
        # heatmap
        heatmap, axis = plt.subplots(figsize=(6,4))
        sns.heatmap(read_counts_df, ax=axis)
        axis.set_title("Read counts")
        heatmap.savefig(output.heatmap, format='svg')
        plt.close()
        # samples barplot
        sample_read_counts = read_counts_df.sum()
        samples = sample_read_counts.plot(kind='bar', title='Read counts per sample', figsize=(6,4))
        plt.savefig(output.barplot_samples, format='svg')
        plt.close()
        # loci barplot
        locus_read_counts = read_counts_df.sum(axis=1)
        loci = locus_read_counts.plot(kind='bar', title='Read counts per locus', figsize=(16,9))
        plt.savefig(output.barplot_loci, format='svg')
        plt.close()

rule count_reads:
    input:
        expand(haps("intermediates_illumina/alignments/{sample}.sorted.bam"), sample=samples)
    output:
        'results/read_counts.csv'
    run:
        read_counts_dict = dict()
        for filename in input:
            sample_id = filename.split('/')[-1].strip('.sorted.bam')
            bam = HTSeq.BAM_Reader(filename)
            read_counts_dict[sample_id] = collections.Counter(aln.iv.chrom for aln in bam if aln.iv)
        read_counts_df = pd.DataFrame(read_counts_dict)
        read_counts_df = read_counts_df.fillna(0)
        read_counts_df = read_counts_df[read_counts_df.columns].astype(int)
        read_counts_df.to_csv(output[0])

rule html_to_pdf:  # TODO: convert plotly html files to pdf. or use seaborn intead of plotly to create SVGs
    input:
        haps_hist=haps("reports/loci_histogram.html"),
        snps_hist=snps("reports/snp_histogram.html")
    output:
        haps_pdf="reports/loci_histogram.pdf",
        snps_pdf="reports/snp_histogram.pdf"
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
        haps_pdf="reports/loci_histogram.pdf",
        haps_mtx=haps("reports/matrix.tsv"),
        snps_pdf="reports/snp_histogram.pdf",
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
