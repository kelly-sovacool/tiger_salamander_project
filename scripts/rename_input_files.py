#!/usr/local/bin/python3
""" Rename raw files to just the sample id & read pair. Strip leading D if followed by numerals.
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
Date: 30 Mar. 2018

Usage:
    ./rename_input_files.py <input_directory>
"""
import os
import re
import sys
input_dir = sys.argv[1]

sample_regex = re.compile('[0-9A-Za-z]*')
leading_D_regex = re.compile('^D[0-9]')
extension_regex = re.compile('\..*')
pair1 = "_R1"
pair2 = "_R2"
for filename in os.listdir(input_dir):
    sample = re.search(sample_regex, filename).group(0)
    if re.search(leading_D_regex, sample):
        sample = sample[1:]  # remove leading D
    extension = re.search(extension_regex, filename).group(0)
    if pair1 in filename:
        pair = pair1
    elif pair2 in filename:
        pair = pair2
    else:
        raise ValueError("neither {} nor {} in {}".format(pair1, pair2, filename))
    new_fn = sample + pair + extension
    if filename != new_fn:
        print('renaming {} to {}'.format(filename, new_fn))
        os.rename(os.path.join(input_dir, filename), os.path.join(input_dir, new_fn))
