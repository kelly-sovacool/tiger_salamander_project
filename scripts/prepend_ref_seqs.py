#!/usr/local/bin/python3
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 20 Feb. 2018
#
# Usage:
#   ./prepend_ref_seqs.py <input_dir> <output_dir> <reference_filename>
import Bio.SeqIO
import collections
import os
import sys


def main():

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    input_dir = input_dir if input_dir.endswith("/") else input_dir + "/"
    output_dir = output_dir if output_dir.endswith("/") else output_dir + "/"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    reference_filename = sys.argv[3]
    input_filenames = os.listdir(input_dir)
    reference_sequences = {seq_record.id: seq_record for seq_record in Bio.SeqIO.parse(reference_filename, 'fasta')}
    for input_filename in input_filenames:
        ref_id = input_filename.split('.')[0]
        if ref_id in reference_sequences:
            ref_seq = reference_sequences[ref_id]
        else:
            raise ValueError("{} not found in reference file.".format(ref_id))
        with open(input_dir + input_filename, 'r') as infile:
            seqs = list(Bio.SeqIO.parse(infile, 'fasta'))
        seqs.insert(0, ref_seq)
        Bio.SeqIO.write(seqs, output_dir + input_filename, 'fasta')

main()
