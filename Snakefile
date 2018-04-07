""" Pipeline for calling SNP haplotypes on Illumina sequences and combining with haplotypes from 454 data. """
# TODO: rule for removing bad loci
# TODO: rule "report" for summarizing matrix
import Bio.SeqIO
import collections
import datetime
import os
import plotly
import time

configfile: "config.yaml"
reference_file = config["reference"]
reference_sequences = {seq_record.id: seq_record for seq_record in Bio.SeqIO.parse(reference_file, 'fasta')}

loci_all=set(reference_sequences.keys())
loci_454 = {fn.strip(".fna") for fn in os.listdir(config["454_haplotypes_dir"])}
loci_bad = set(config["bad_loci"])
loci_bad.update(loci_all - loci_454)
loci = loci_all - loci_bad

pair1, pair2 = "R1", "R2"
hap_numbers = ["1", "2"]

include: "rules/haplotype_illumina_data.smk"
include: "rules/align_454_haplotypes.smk"

rule all:
    input:
        config["matrix_filename"],
        expand("snp_sites/{locus}.fna", locus=loci)

rule combine_454_and_illumina:
    input:
        haps_454="intermediates_454/haps_aligned_filtered/{locus}.fna",
        haps_illumina="intermediates_illumina/haplotypes/{locus}.fna"
    output:
        protected("haplotypes/{locus}.fna")
    shell:
        "cat {input} > {output}"

rule remove_bad_loci:
    input:
        expand("haplotypes/{locus}.fna", locus=loci)
    output:
        "logs/removed_bad_loci.txt"
    run:
        removed_loci = list()
        for fn in input:
            locus = fn.split('/')[1].split('.')[0]
            if locus in bad_loci:
                os.remove(fn)
                removed_loci.append(locus)
        with open(output[0], 'w') as file:
            file.write('\n'.join(removed_loci))

rule snp_sites:
    input:
        "haplotypes/{locus}.fna"
    output:
        "snp_sites/{locus}.fna"
    log:
        "logs/snp_sites/{locus}.log"
    benchmark:
        "benchmarks/snp_sites/{locus}.txt"
    shell:
        "snp-sites -o {output} {input} 2> {log}"

rule build_matrix:
    input:
        expand("haplotypes/{locus}.fna", locus=loci)
    output:
        matrix=config["matrix_filename"],
        histogram=config["hist_filename"]
    benchmark:
        "benchmarks/build_matrix.txt"
    run:
        individuals_to_seq_runs = collections.defaultdict(set)
        for working_dir, subdirs, files in os.walk(os.path.expanduser(config['illumina_fastq_dir'])):
            for filename in files:
                if filename[:4] != "Icon":
                    indiv_id = filename.split('_')[0]
                    base, subdir = os.path.split(working_dir)
                    individuals_to_seq_runs[indiv_id].add(subdir)
        for filename in os.listdir(config["454_haplotypes_dir"]):
            with open(os.path.join(config["454_haplotypes_dir"], filename), 'r') as infile:
                for seq_record in Bio.SeqIO.parse(infile, 'fasta'):
                    individuals_to_seq_runs[seq_record.id].add("454")
        matrix = "individual_id\tsequencing_runs"
        input = sorted(input)
        loci = [filename.split('/')[-1].split('.fna')[0] for filename in input]
        individuals_to_loci = collections.defaultdict(set)
        for locus_filename in input:
            locus_id = locus_filename.split('/')[-1].split('.fna')[0]
            matrix += "\t" + locus_id
            with open(locus_filename, 'r') as file:
                for line in file:
                    if line[0] == '>':
                        indiv_id = line[1:].split('_')[0]
                        individuals_to_loci[indiv_id].add(locus_id)
        for indiv_id in sorted(individuals_to_loci):
            matrix += '\n' + indiv_id + '\t' + ', '.join(sorted(individuals_to_seq_runs[indiv_id]))
            for locus_id in loci:
                matrix += '\t1' if locus_id in individuals_to_loci[indiv_id] else '\t0'
        with open(output.matrix, 'w') as file:
            file.write(matrix)
        plotly.offline.plot(plotly.graph_objs.Figure(data=[plotly.graph_objs.Histogram(x=len(individuals_in_loci[indiv_id]), name=indiv_id, opacity=0.75, autobinx=True) for indiv_id in sorted(individuals_in_loci)],
                            layout=plotly.graph_objs.Layout(barmode='overlay', title=('Individuals in Loci'), xaxis=dict(title='Number of Loci'), yaxis=dict(title='Number of Individuals'))),
                            filename=output.histogram, auto_open=True)

rule snp_subsample:
    input:
        expand("snp_sites/{locus}.fna", locus=loci)
    output:
        expand("snp_subsamples/{subsample_base}{num}.{{ext}}", subsample_base=config["snp_subsample"]["output_filename_base"], num=range(config["snp_subsample"]["num_subsamples"])),
        "snp_subsample/all-snps-all-loci.{ext}"
    params:
        fn_base = config["snp_subsample"]["output_filename_base"],
        num_subsamples = config["snp_subsample"]["num_subsamples"],
        output_format = config["snp_subsample"]["output_format"]
    shell:
        "scripts/snp_subsample.py snp_sites/ snp_subsamples/{params.fn_base} --all-snps-all-loci "
        "--output-format {params.output_format} --num_subsamples {params.num_subsamples}"

rule write_strx_scripts:
    input:
        "snp_subsamples/{subsample_filename}.{ext}"
    output:
        expand("structure/strx_{{subsample_filename}}_{k}.sh", k=range(config["structure"]["k_min"], config["structure"]["k_max"] + 1))
    params:
        strx_path = config["structure"]["path"],
        num_runs = config["structure"]["num_runs"],
        k_min = config["structure"]["k_min"],
        k_max = config["structure"]["k_max"],
        email = config["structure"]["email"]
    run:
        N = 0  # individuals
        L = 0  # loci
        with open(input[0], 'r') as file:
            is_first = True
            for line in file:
                N += 1
                if is_first:
                    is_first = False
                    L = len(line.split()) - 1
        N = N / 2
        for k in range(params.k_min, params.k_max + 1):
            script = "#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user={email}\nSBATCH_NODELIST: $SBATCH_NODELIST\n\n".format(email=params.email)
            for run in range(1, params.num_runs + 1):
                script += "{strx_path} -i {subsample} -K {k} -L {L} -N {N} -o {{subsample_filename}}_k{k}_run{run}.out &\nsleep 10\n\n".format(strx_path=params.strx_path, subsample=input[0], k=k, L=L, N=N, run=run)
            script += 'wait'
            with open('structure/strx_{{subsample}}_k{k}.sh'.format(k=k), 'w') as script_file:
                script_file.write(script)

rule report:
    input:  # TODO: report snps per locus histogram
        snps=expand("snp_sites/{locus}.fna", locus=loci),
        matrix=config["matrix_filename"],
        hist=config["hist_filename"],
        removed_loci_file="logs/removed_bad_loci.txt"
    output:  # TODO: histogram showing number of individuals per locus
        "reports/report.html"
    run:
        from snakemake.utils import report
        with open(input.matrix) as file:
            matrix = file.readlines()
        num_individuals = len(matrix) - 1
        num_loci = len(loci)
        with open(input.removed_bad_loci, 'r') as file:
            removed_loci = file.readlines()
        report("""
        SNP calling & haplotyping pipeline for the Tiger Salamander Project
        ===================================

        {num_individuals} individuals were haplotyped for {num_loci} loci (see Table T1_ and Graph G1_).
        Loci removed: {removed_loci}
        """, output[0], T1=input.matrix, G1=input.hist)
