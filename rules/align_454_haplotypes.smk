""" Pipeline for processing 454 haplotypes before combining with Illumina data """
import Bio.SeqIO
import collections
import os

input_dir = config["454_haplotypes_dir"]
wildcards = glob_wildcards("{haplotypes_dir}/{{locus}}.fna".format(haplotypes_dir=input_dir))
loci = wildcards.locus

rule insert_references:
    input:
        "{haplotypes_dir}/{{locus}}.fna".format(haplotypes_dir=input_dir)
    output:
        "intermediates_454/haps_combined_with_refs/{locus}.fna"
    params:
        locus="{locus}"
    run:
        with open(input[0], 'r') as infile:
            sequences = list(Bio.SeqIO.parse(infile, 'fasta'))
        sequences.insert(0, reference_sequences[params.locus])
        with open(output[0], 'w') as outfile:
            Bio.SeqIO.write(sequences, outfile, 'fasta')

rule align:
    input:
        "intermediates_454/haps_combined_with_refs/{locus}.fna"
    output:
        "intermediates_454/haps_aligned_with_refs/{locus}.fna"
    shell:
        "mafft --auto --quiet {input} > {output}"

rule filter:
    input:
        "intermediates_454/haps_aligned_with_refs/{locus}.fna"
    output:
        protected("intermediates_454/haps_aligned_filtered/{locus}.fna")
    params:
        locus="{locus}"
    run:
        nucleotides = "agtcAGTC"
        output_sequences = list()
        input_sequences = {seq_record.id: seq_record for seq_record in Bio.SeqIO.parse(input[0], 'fasta')}
        locus_record = input_sequences.pop(params.locus)  # remove reference sequence
        insertion_positions = {pos for pos, nuc in enumerate(locus_record.seq) if nuc == '-'}
        for seq_id, seq_record in input_sequences.items():
            # replace ambiguous nucleotides with dashes
            new_seq = "".join(nuc if nuc in nucleotides else '-' for pos, nuc in enumerate(seq_record.seq) if pos not in insertion_positions)
            if ''.join(set(new_seq)) != '-': # throw out sequences that are entirely dashes
                seq_record.seq = Bio.Seq.Seq(new_seq)
                output_sequences.append(seq_record)
            else:
                with open(log[0], 'w') as logfile:
                    logfile.write("{} from {} has no nucleotides".format(new_seq.id, input[0]))
        Bio.SeqIO.write(output_sequences, output[0], 'fasta')
