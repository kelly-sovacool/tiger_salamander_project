#!/usr/local/bin/python3
""" Remove an individual from all fasta files in a directory.
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
31 Jan. 2018

Usage:
    ./remove_individual.py <individual_id> <directory>

Options:
    -h --help   Display this help message.
"""
import Bio.SeqIO
import docopt
import os
import shutil


def main(args):
    directory = args['<directory>']
    for thing in os.listdir(directory):
        if not thing.endswith("/"):
            filename = directory + thing
            backup = filename + '.bak'
            shutil.copy2(filename, backup)
            with open(filename, 'r') as input_file:
                records = [rec for rec in Bio.SeqIO.parse(input_file, 'fasta') if rec.id.split('_')[0] != args['<individual_id>']]
            Bio.SeqIO.write(records, filename, 'fasta')


def check_directory(dir_name):
        if not os.path.isdir(dir_name):
            raise FileNotFoundError("Directory {} does not exist.".format(dir_name))
        if not dir_name.endswith('/'):
            dir_name = dir_name + '/'
        return dir_name


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    for arg in arguments:
        if arguments[arg] and "dir" in arg:
            arguments[arg] = check_directory(arguments[arg])
    main(arguments)

