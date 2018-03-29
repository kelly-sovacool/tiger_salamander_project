#!/usr/local/bin/python3
""" Combine fasta files from multiple sequencing runs to make unpaired fasta files, one per individual
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
Date: 29 Mar. 2018

Usage:
    ./combine_seq_runs.py <seq_runs_dir> <output_dir>

Options:
    -h --help   display this super helpful message
"""
import Bio.SeqIO
import collections
import docopt
import os


def main(args):
    individual_ids = set()
    if not os.path.isdir(args['<output_dir>']):
        os.mkdir(args['<output_dir>'])
    for subdir in os.listdir(args['<seq_runs_dir>']):
        if os.path.isdir(args['<seq_runs_dir>'] + subdir):
            subdir = check_directory(subdir)
            for filename in os.listdir(args['<seq_runs_dir>'] + subdir):
                if os.path.isfile(args['<seq_runs_dir>'] + subdir + filename):
                    indiv_id = filename.split('_')[0]
                    individual_ids.add(indiv_id)
                    with open(args['<seq_runs_dir>'] + subdir + filename, 'r') as input_file:
                        with open(args['<output_dir>'] + indiv_id + '.tmp', 'a') as output_file:
                            for line in input_file:
                                output_file.write(line)
    for indiv_id in individual_ids:
        with open(args['<output_dir>'] + indiv_id + '.tmp', 'r') as input_file:
            with open(args['<output_dir>'] + indiv_id + '.fna', 'w') as output_file:
                count = 0
                for line in input_file:
                    if line[0] == '>':
                        count += 1
                        line = '>' + indiv_id + '_' + count
                    output_file.write(line)
        os.remove(args['<output_dir>'] + indiv_id + '.tmp')


def check_directory(dir_name):
        if not dir_name.endswith('/'):
            dir_name = dir_name + '/'
        return dir_name

if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    for arg in arguments:
        arguments[arg] = check_directory(arguments[arg])
    main(arguments)
