""" Pipeline for processing Illumina MiSeq amplicon sequence data and calling SNP haplotypes. """

import Bio.SeqIO
import collections
import os
import re
import shutil

fastq_dir = config["illumina_fastq_dir"]
adapter_file = config['trimmomatic_adapter']
loci=set(reference_sequences.keys())

def get_sample_ids(input_dir, r1, r2):
    """ Renames files to strip the leading D and keep only the id & read pair.
    Warning: do not use this if any individual/sample was sequenced on multiple lanes.
    Return: dict mapping sample ids to list of raw fastq filepaths (e.g. for individuals in multiple sequencing runs)
    """  # TODO: support samples on multiple lanes
    samples_to_seq_runs = collections.defaultdict(dict)
    sample_regex = re.compile('[0-9A-Za-z]*')
    leading_D_regex = re.compile('^D[0-9]')
    extension_regex = re.compile('\..*')
    for working_dir, subdir, files in os.walk(os.path.expanduser(input_dir)):
        for filename in files:
            if filename[0] != '.' and filename[:4] != "Icon" :  # don't work on hidden files
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

samples_to_filepaths = get_sample_ids(fastq_dir, pair1, pair2)  # sample_id: [seq_run1/id.fastq.gz, seq_run2/id.fastq.gz, ...]
wildcards = glob_wildcards(os.path.join(fastq_dir, "{seq_run}/{sample}_{pair}.fastq.gz"))
samples = set(wildcards.sample)

rule haplotype_illumina:
    input:
        expand("intermediates_illumina/individual_haps/{sample}.hap{h}.fna", sample=samples, h=hap_numbers),
        expand("intermediates_illumina/haplotypes/{locus}.fna", locus=loci)

rule fastq_merge:
    input:
        lambda wildcards: samples_to_filepaths[wildcards.sample][wildcards.pair]
    output:
        "intermediates_illumina/merged_fastq/{sample}_{pair}.fastq.gz"
    log:
        "logs/fastq_merge/{sample}_{pair}.log"
    run:
        if len(input) == 1:
            shutil.copy2(input[0], output[0])
        else:
            with open(log[0], 'w') as logfile:
                logfile.write("Merging {}".format(str(input)))
            with open(output[0], 'wb') as out_file:
                for filename in input:
                    with open(filename, 'rb') as in_file:
                        shutil.copyfileobj(in_file, out_file)

rule fastq_filter:
    input:
        R1="intermediates_illumina/merged_fastq/{{sample}}_{pair}.fastq.gz".format(pair=pair1),
        R2="intermediates_illumina/merged_fastq/{{sample}}_{pair}.fastq.gz".format(pair=pair2)
    output:
        R1_paired="intermediates_illumina/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair1),
        R1_single="intermediates_illumina/filtered_fq/{{sample}}_{pair}_single.fastq.gz".format(pair=pair1),
        R2_paired="intermediates_illumina/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair2),
        R2_single="intermediates_illumina/filtered_fq/{{sample}}_{pair}_single.fastq.gz".format(pair=pair2)
    log:
        "logs/fastq_filter/{sample}.log"
    shell:
        "trimmomatic PE -phred33 {input.R1} {input.R2} {output.R1_paired} {output.R1_single} {output.R2_paired} {output.R2_single} "
        "ILLUMINACLIP:{adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &> {log}"

rule bwa_index:
    input:
        reference_file
    output:
        reference_file + '.bwt'
    log:
        "logs/bwa_index_{ref}.log".format(ref=reference_file)
    shell:
        "bwa index {input} &> {log}"

rule bwa_map:
    input:
        R1="intermediates_illumina/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair1),
        R2="intermediates_illumina/filtered_fq/{{sample}}_{pair}.fastq.gz".format(pair=pair2),
        ref=reference_file + '.bwt'
    output:
        temp("intermediates_illumina/alignments/{sample}.mapped.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_map/{sample}.log"
    benchmark:
        "benchmarks/bwa_map/{sample}.txt"
    shell:
        "bwa mem -R '{params.rg}' {reference_file} {input.R1} {input.R2} | "
        "samtools view -Sb - > {output} 2> {log}"

rule samtools_sort:
    input:
        "intermediates_illumina/alignments/{sample}.mapped.bam"
    output:
        "intermediates_illumina/alignments/{sample}.sorted.bam"
    benchmark:
        "benchmarks/samtools_sort/{sample}.txt"
    shell:
        "samtools sort -T intermediates_illumina/alignments/{wildcards.sample} -O bam {input} > {output}"

rule samtools_index:
    input:
        "intermediates_illumina/alignments/{sample}.sorted.bam"
    output:
        "intermediates_illumina/alignments/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_faidx:
    input:
        reference_file
    output:
        reference_file + '.fai'
    shell:
        "samtools faidx {input}"

rule call_variants:
    input:
        bam="intermediates_illumina/alignments/{sample}.sorted.bam",
        fai=reference_file + '.fai'
    output:
        "intermediates_illumina/variants/calls/{sample}.vcf"
    log:
        "logs/call_variants/{sample}.log"
    benchmark:
        "benchmarks/call_variants/{sample}.txt"
    shell:
        "freebayes --min-mapping-quality 1 -f {reference_file} {input.bam} > {output} 2> {log}"

rule filter_variants:
    input:
        "intermediates_illumina/variants/calls/{sample}.vcf"
    output:
        "intermediates_illumina/variants/filtered/{sample}.vcf"
    log:
        "logs/filter_variants/{sample}.log"
    benchmark:
        "benchmarks/filter_variants/{sample}.txt"
    shell:
        "vcftools --remove-indels --vcf {input} --max-missing 0.5 --minQ 20 --minDP 3 --recode --recode-INFO-all -c > {output} 2> {log}"

rule phase_variants:
    input:
        vcf="intermediates_illumina/variants/filtered/{sample}.vcf",
        bam="intermediates_illumina/alignments/{sample}.sorted.bam",
        bai="intermediates_illumina/alignments/{sample}.sorted.bam.bai"
    output:
        "intermediates_illumina/variants/phased/{sample}.vcf"
    log:
        "logs/phase_variants/{sample}.log"
    benchmark:
        "benchmarks/phase_variants/{sample}.txt"
    shell:
        "whatshap phase --reference {reference_file} -o {output} {input.vcf} {input.bam} &> {log}"

rule bgzip:
    input:
        "intermediates_illumina/variants/phased/{sample}.vcf"
    output:
        "intermediates_illumina/variants/phased/{sample}.vcf.gz"
    shell:
        "bgzip {input}"

rule tabix:
    input:
        "intermediates_illumina/variants/phased/{sample}.vcf.gz"
    output:
        "intermediates_illumina/variants/phased/{sample}.vcf.gz.tbi"
    shell:
        "tabix {input}"

rule consensus:
    input:
        vcf="intermediates_illumina/variants/phased/{sample}.vcf.gz",
        tbi="intermediates_illumina/variants/phased/{sample}.vcf.gz.tbi",
        bam="intermediates_illumina/alignments/{sample}.sorted.bam"
    output:
        "intermediates_illumina/individual_haps/{sample}.hap{h}.fna"
    params:
        sample="{sample}",
        hap_number="{h}"
    shell:
        "bcftools consensus -i -s {params.sample} -H {params.hap_number} -f {reference_file} {input.vcf} > {output}"

rule combine_haps:
    input:
        expand("intermediates_illumina/individual_haps/{sample}.hap{h}.fna", sample=samples, h=hap_numbers)
    output:
        expand("intermediates_illumina/haps_combined_with_refs/{locus}.fna", locus=loci)
    benchmark:
        "benchmarks/combine_haps.txt"
    run:
        haps_in_loci=collections.defaultdict(dict)
        for sample in samples:
            for filename, hap_number in (("intermediates_illumina/individual_haps/{}.hap1.fna".format(sample), '_1'), ("intermediates_illumina/individual_haps/{}.hap2.fna".format(sample), '_2')):
                for seq_record in Bio.SeqIO.parse(filename, 'fasta'):
                    locus_id = seq_record.id
                    seq_record.id = sample + hap_number
                    haps_in_loci[locus_id][seq_record.id] = seq_record
        for locus_id in haps_in_loci:
            sequences = [haps_in_loci[locus_id][seq_id] for seq_id in sorted(haps_in_loci[locus_id].keys())]
            sequences.insert(0, reference_sequences[locus_id])
            with open("intermediates_illumina/haps_combined_with_refs/" + locus_id + '.fna', 'w') as file:
                Bio.SeqIO.write(sequences, file, 'fasta')

rule align_haps:
    input:
        "intermediates_illumina/haps_combined_with_refs/{locus}.fna"
    output:
        "intermediates_illumina/haps_aligned_with_refs/{locus}.fna"
    wildcard_constraints:
        locus="\w+"
    benchmark:
        "benchmarks/align_haps/{locus}.txt"
    shell:
        "mafft --auto --quiet {input} > {output}"

rule filter_haps:
    input:
        "intermediates_illumina/haps_aligned_with_refs/{locus}.fna"
    output:
        "intermediates_illumina/haplotypes/{locus}.fna"
    wildcard_constraints:
        locus="\w+"
    log:
        "logs/output_haps/{locus}.log"
    benchmark:
        "benchmarks/filter_haps/{locus}.txt"
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
                    with open(log[0], 'w') as logfile:
                        logfile.write("{} from {} has no nucleotides".format(new_seq.id, input[0]))
        Bio.SeqIO.write(sequences, output[0], 'fasta')
