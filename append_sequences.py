#!/usr/local/bin/python3
""" An ad hoc script to add sequences from our 454 run to our illumina runs, 
    throwing out 454 individuals that are already in the illumina datasets.
    Also sorts the output by sequence ids.

Author: Kelly Sovacool
Email: kellysovacool@uky.edu
Date: 18 Feb. 2018

Usage:
    ./append_sequences.py <input_path> <output_path>

Options:
    -h --help   Display this dubiously helpful message
"""
import Bio.SeqIO
import docopt
import os


def main(args):
    input_path = args['<input_path>']
    output_path = args['<output_path>']
    if os.path.isdir(input_path):
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        input_path = input_path if input_path.endswith("/") else input_path + "/"
        output_path = output_path if output_path.endswith("/") else output_path + "/"
        for input_filename in os.listdir(input_path):
            append_sequences(input_path + input_filename, output_path + input_filename)
    elif os.path.isfile(input_path):
        append_sequences(input_path, output_path)
    else:
        print("Cannot do conversion. Input {} and output {} must either both be files or both be directories.".format(input_path, output_path))


def append_sequences(input_filename, output_filename):
    with open(output_filename, 'r') as out_file:
        all_seqs = list(Bio.SeqIO.parse(out_file, 'fasta'))
        ids = {seq.id for seq in all_seqs}
    with open(input_filename, 'r') as in_file:
        for new_seq in Bio.SeqIO.parse(in_file, 'fasta'):
            if new_seq.id not in ids:
                all_seqs.append(new_seq)
            else:
                print(new_seq.id, "already in", output_filename)
    sorted_seqs= sorted(all_seqs, key=lambda seq: seq.id)
    Bio.SeqIO.write(sorted_seqs, output_filename, 'fasta')

if __name__ == "__main__":
    main(docopt.docopt(__doc__))

