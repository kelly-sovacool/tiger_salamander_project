""" Pipeline for processing data from amplicon sequencing and calling haplotypes.
Author: Kelly Sovacool
Email: kellysovacool@gmail.com
Date: 30 Mar. 2018
"""
# TODO: support for multiple sequencing runs, but optional. use samtools merge
# TODO: use log files
# TODO: use benchmarking

import Bio.SeqIO
import collections
import os
import re
import shutil

configfile: "config.yaml"

fastq_dir = config["fastq_dir"]
#seq_runs = config["seq_runs"]
haps_dir=config["haps_dir"]
reference_file = config["reference"]
reference_sequences = {seq_record.id: seq_record for seq_record in Bio.SeqIO.parse(reference_file, 'fasta')}
loci=set(reference_sequences.keys())
jar_file = config["trimmomatic_jar"]
adapter_file = config['trimmomatic_adapter']
pair1, pair2 = "R1", "R2" # set of paired end specifiers

def get_sample_ids(input_dir, r1, r2):
    samples_to_seq_runs = collections.defaultdict(dict)
    sample_regex = re.compile('[0-9A-Za-z]*')
    leading_D_regex = re.compile('^D[0-9]')
    extension_regex = re.compile('\..*')
    for working_dir, subdir, files in os.walk(os.path.expanduser(input_dir)):
        for filename in files:
            if filename[0] != '.':  # don't work on hidden files
                sample = re.search(sample_regex, filename).group(0)
                if re.search(leading_D_regex, sample):
                    sample = sample[1:]  # remove leading D
                extension = re.search(extension_regex, filename).group(0)
                if "_" + r1 in filename:
                    pair = r1
                elif "_" + r2 in filename:
                    pair = r2
                else:
                    raise ValueError("neither {} nor {} in {}".format(r1, r2, filename))
                old_filepath = os.path.join(working_dir, filename)
                new_filepath =  os.path.join(working_dir, sample + "_" + pair + extension)
                if old_filepath != new_filepath:
                    os.rename(old_filepath, new_filepath)
                    current_filepath = new_filepath
                else:
                    current_filepath = old_filepath
                if pair not in samples_to_seq_runs[sample]:
                    samples_to_seq_runs[sample][pair] = list()
                samples_to_seq_runs[sample][pair].append(current_filepath)
    return samples_to_seq_runs

samples_to_filepaths = get_sample_ids(fastq_dir, pair1, pair2)  # sample_id: {seq_run1/id.fastq.gz, seq_run2/id.fastq.gz, ...}
wildcards = glob_wildcards(os.path.join(fastq_dir, "{seq_run}/{sample}_{pair}.fastq.gz"))
samples = set(wildcards.sample)

rule all:
    input:
        "matrix.tsv"

rule fastq_merge:
    input:
        lambda wildcards: samples_to_filepaths[wildcards.sample][wildcards.pair]
    output:
        "intermediates/merged_fastq/{sample}_{pair}.fastq.gz"
    run:
        if len(input) == 1:
            shutil.copy2(input[0], output[0])
        else:
            with open(output[0], 'wb') as out_file:
                for filename in input:
                    with open(filename, 'rb') as in_file:
                        shutil.copyfileobj(in_file, out_file)

rule fastq_filter:
    input:
        R1="intermediates/merged_fastq/{{sample}}_{pair}.fastq.gz".format(pair=pair1),
        R2="intermediates/merged_fastq/{{sample}}_{pair}.fastq.gz".format(pair=pair2)
    output:
        R1_paired="intermediates/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair1),
        R1_single="intermediates/filtered_fq/{{sample}}_{pair}_single.fastq.gz".format(pair=pair1),
        R2_paired="intermediates/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair2),
        R2_single="intermediates/filtered_fq/{{sample}}_{pair}_single.fastq.gz".format(pair=pair2)
    shell:
        "java -jar {jar_file} PE -phred33 {input.R1} {input.R2} {output.R1_paired} {output.R1_single} {output.R2_paired} {output.R2_single} "
        "ILLUMINACLIP:{adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule bwa_index:
    input:
        reference_file
    output:
        reference_file + '.bwt'
    shell:
        "bwa index {input}"

rule bwa_map:
    input:
        R1="intermediates/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair1),
        R2="intermediates/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair2),
        ref=reference_file + '.bwt'
    output:
        "intermediates/alignments/{sample}.mapped.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "bwa mem {reference_file} {input.R1} {input.R2} | "
        "samtools view -Sb - > {output} 2> {log}"

rule samtools_sort:
    input:
        "intermediates/alignments/{sample}.mapped.bam"
    output:
        protected("intermediates/alignments/{sample}.sorted.bam")
    shell:
        "samtools sort -T intermediates/alignments/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_faidx:
    input:
        reference_file
    output:
        reference_file + '.fai'
    shell:
        "samtools faidx {input}"

rule freebayes:
    input:
        bam="intermediates/alignments/{sample}.sorted.bam",
        ref=reference_file + '.fai'
    output:
        "intermediates/variants/{sample}.raw.vcf"
    shell:
        "freebayes --min-mapping-quality 1 -f {reference_file} {input.bam} > {output}"

rule bcftools_filter:
    input:
        "intermediates/variants/{sample}.raw.vcf"
    output:
        "intermediates/variants/{sample}.filtered.vcf"
    shell:
        "bcftools filter -e 'QUAL < 20' {input} | bcftools filter -e 'DP < 5' -o {output}"

rule whatshap:
    input:
        vcf="intermediates/variants/{sample}.filtered.vcf",
        bam="intermediates/alignments/{sample}.sorted.bam"
    output:
        "intermediates/variants/{sample}.phased.vcf"
    shell:
        "whatshap phase --ignore-read-groups --reference {reference_file} -o {output} {input.vcf} {input.bam}"

rule bgzip:
    input:
        "intermediates/variants/{sample}.phased.vcf"
    output:
        "intermediates/variants/{sample}.phased.vcf.gz"
    shell:
        "bgzip {input}"

rule tabix:
    input:
        "intermediates/variants/{sample}.phased.vcf.gz"
    output:
        "intermediates/variants/{sample}.phased.vcf.gz.tbi"
    shell:
        "tabix {input}"

rule consensus1:
    input:
        vcf="intermediates/variants/{sample}.phased.vcf.gz",
        tbi="intermediates/variants/{sample}.phased.vcf.gz.tbi"
    output:
        "intermediates/individual_haps/{sample}.hap1.fna"
    shell:
        "bcftools consensus -H 1 -f {reference_file} {input.vcf} > {output}"

rule consensus2:
    input:
        vcf="intermediates/variants/{sample}.phased.vcf.gz",
        tbi="intermediates/variants/{sample}.phased.vcf.gz.tbi"
    output:
        "intermediates/individual_haps/{sample}.hap2.fna"
    shell:
        "bcftools consensus -H 2 -f {reference_file} {input.vcf} > {output}"

rule combine_haps:
    input:
        expand("intermediates/individual_haps/{sample}.hap{h}.fna", sample=samples, h={1,2})
    output:
        expand("intermediates/combined_haps_with_refs/{locus}.haps.fna", locus=loci)
    run:
        haps_in_loci=collections.defaultdict(dict)
        for sample in samples:
            for filename, hap_number in (("intermediates/individual_haps/{}.hap1.fna".format(sample), '_1'), ("intermediates/individual_haps/{}.hap2.fna".format(sample), '_2')):
                for seq_record in Bio.SeqIO.parse(filename, 'fasta'):
                    locus_id = seq_record.id
                    seq_record.id = sample + hap_number
                    haps_in_loci[locus_id][seq_record.id] = seq_record
        for locus_id in haps_in_loci:
            sequences = [haps_in_loci[locus_id][seq_id] for seq_id in sorted(haps_in_loci[locus_id].keys())]
            sequences.insert(0, reference_sequences[locus_id])
            with open("intermediates/combined_haps_with_refs/" + locus_id + '.haps.fna', 'w') as file:
                Bio.SeqIO.write(sequences, file, 'fasta')

rule align_haps:
    input:
        "intermediates/combined_haps_with_refs/{locus}.haps.fna"
    output:
        "intermediates/aligned_hap_with_refs/{locus}.haps.fna"
    wildcard_constraints:
        locus="\w+"
    shell:
        "mafft --retree 1 --maxiterate 0 --quiet {input} > {output}"

rule output_haps:
    input:
        "intermediates/aligned_hap_with_refs/{locus}.haps.fna"
    output:
        "{haps_dir}/{locus}.fna"
    wildcard_constraints:
        locus="\w+"
    run:
        nucleotides = "agtcAGTC"
        sequences = list()
        for seq_record in Bio.SeqIO.parse(input[0], 'fasta'):
            # remove reference sequences
            if seq_record.id not in loci:
                # replace ambiguous nucleotides with dashes
                new_seq = "".join(nuc if nuc in nucleotides else '-' for nuc in seq_record.seq)
                if ''.join(set(new_seq)) != '-': # throw out sequences that are entirely dashes
                    seq_record.seq = Bio.Seq.Seq(new_seq)
                    sequences.append(seq_record)
                else:
                    print("{} from {} has no nucleotides".format(new_seq.id, input[0]))
        Bio.SeqIO.write(sequences, output[0], 'fasta')

rule build_matrix:
    input:
        expand("{haps_dir}/{locus}.fna", locus=loci, haps_dir=haps_dir)
    output:
        "matrix.tsv"
    run:
        matrix = "individual_id"
        loci = list(filename.split('/')[-1].split('.fna')[0] for filename in sorted(input))
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
            matrix += '\n' + indiv_id
            for locus_id in loci:
                matrix += '\t1' if locus_id in individuals_to_loci[indiv_id] else '\t0'
        with open(output[0], 'w') as file:
            file.write(matrix)
