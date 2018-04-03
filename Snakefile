""" Pipeline for processing amplicon sequence data and calling SNP haplotypes.
Author: Kelly Sovacool
Date: 30 Mar. 2018
"""
# TODO: remove bad locus

import Bio.SeqIO
import collections
import os

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
        config["matrix_filename"]
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
        with open(output[0], 'w') as file:
            file.write(matrix)
