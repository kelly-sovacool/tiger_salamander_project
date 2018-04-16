#!/usr/local/bin/python3

import Bio.SeqIO
import sys

input_filename = sys.argv[1]

with open(input_filename, 'r') as infile:
    first = True
    prev = None
    count = 0
    for record in Bio.SeqIO.parse(infile, 'fasta'):
        count += 1
        if first:
            prev = record
            first = False
        else:
            first = True
            if prev.id.split('_')[0] != record.id.split('_')[0]:
                print('widowed', prev.id)
                break

print("total sequences:", count)
