#!/usr/local/bin/python3
"""
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
Date: 27 Mar. 2018

Usage: build_ids_in_loci_mtx.py <loci_dir> <seq_runs_dir> <output_filename>

Options:
    -h --help   display this helpful message
"""
import collections
import docopt
import os
import pprint


def main(args):
    loci_names = [name.replace("_ids.txt", '') for name in sorted(os.listdir(args['<loci_dir>']))]
    matrix = "individual_id\tsequencing_run"
    for locus_name in loci_names:
        matrix += "\t" + locus_name
    matrix += '\n'
    indiv_ids_in_loci = collections.defaultdict(set)
    for locus_filename in os.listdir(args['<loci_dir>']):
        locus_name = locus_filename.replace("_ids.txt", '')
        with open(args['<loci_dir>'] + locus_filename, 'r') as file:
            indiv_ids_in_loci[locus_name].update(indiv_id for indiv_id in file)
    for seq_run_filename in os.listdir(args['<seq_runs_dir>']):
        seq_run_name = str(seq_run_filename).strip().replace("_indiv_ids.txt", '')
        with open(args['<seq_runs_dir>'] + seq_run_filename, 'r') as file:
            for indiv_id in file:
                output_line = str(indiv_id).strip() + "\t" + seq_run_name
                for locus in loci_names:
                    if indiv_id in indiv_ids_in_loci[locus]:
                        output_line += "\t1"
                    else:
                        output_line += "\t0"
                matrix += output_line + '\n'
    with open(args['<output_filename>'], 'w') as file:
        file.write(matrix)


if __name__ == "__main__":
    main(docopt.docopt(__doc__))