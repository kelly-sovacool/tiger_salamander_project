
import Bio.SeqIO
import collections
import os

outdir="intermediates/individuals/"
haplotypes_dir="data/"

samples = {seq_record.id.split('_')[0] for locus_filename in os.listdir(haplotypes_dir) for seq_record in Bio.SeqIO.parse(os.path.join(haplotypes_dir, locus_filename), 'fasta')}
loci = {fn.strip('.fna') for fn in os.listdir(haplotypes_dir)}

if samples != {fn.strip('.fna') for fn in os.listdir(outdir)}:
    individuals = collections.defaultdict(list)
    for locus_filename in os.listdir(haplotypes_dir):
        for seq_record in Bio.SeqIO.parse(os.path.join(haplotypes_dir, locus_filename), 'fasta'):
            sample = seq_record.id.split('_')[0]
            individuals[sample].append(seq_record)
    for sample, sequences in individuals.items():
        Bio.SeqIO.write(sequences, "intermediates/individuals/{}.fna".format(sample), 'fasta')
