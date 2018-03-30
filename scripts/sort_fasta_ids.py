#!/usr/local/bin/python3
"""
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
07 Feb. 2018

Usage:
    ./sort_fasta_ids.py input_dir output_dir
"""

import Bio.SeqIO
import os
import sys

def main():
    input_dir = check_directory(sys.argv[1])
    output_dir = check_directory(sys.argv[2])
    for infilename in os.listdir(input_dir):
        with open(input_dir + infilename, 'r') as infile:
            sequences = sorted(Bio.SeqIO.parse(infile, 'fasta'), key=lambda seq: seq.id)
        Bio.SeqIO.write(sequences, output_dir + infilename, 'fasta')

def check_directory(dir_name):
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        if not dir_name.endswith('/'):
            dir_name = dir_name + '/'
        return dir_name

main()
