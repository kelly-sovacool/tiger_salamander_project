#!/usr/local/bin/python3
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 28 Sept. 2017
#
# Report duplicate lines and write non-duplicates to a new file.
# Usage:
#   ./report_duplicate_lines.py <input_filename> <output_filename>

import sys

def main():
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    seq_ids = set()
    with open(input_filename, 'r') as infile:
        for line in infile:
            seq_id = line.strip()
            if seq_id in seq_ids:
                print("duplicate:", seq_id)
            seq_ids.add(seq_id)
    with open(output_filename, "w") as outfile:
        for seq_id in sorted(seq_ids):
            outfile.write(seq_id + "\n")

main()
