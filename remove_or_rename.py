#!/usr/local/bin/python3

import os
import subprocess

filenames = set(os.listdir('./'))
for fn in filenames:
    if fn.endswith('.fasta'):
        if os.path.isfile(fn.replace('.fasta', '.fna')):
            os.remove(fn)
        else:
            subprocess.call(['fastx_renamer', '-n', 'COUNT', '-i', fn, '-o', fn.replace('.fasta', '.fna')])
            os.remove(fn)
