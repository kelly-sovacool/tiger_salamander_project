#!/usr/local/bin
"""
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
28 Sept. 2017

Usage:
    count_blast_records_and_hits.py <blast_result_xml>
"""

import Bio.Blast.NCBIXML
import docopt


def main(args):
    count_records = 0
    count_alignments = 0
    count_hits = 0
    with open(args['<blast_result_xml>'], 'r') as xml_file:
        for record in Bio.Blast.NCBIXML.parse(xml_file):
            count_records += 1
            for alignment in record.alignments:
                count_alignments += 1
                for hit in alignment.hsps:
                    count_hits += 1
    print('records:', count_records)
    print('alignments:', count_alignments)
    print('count_hits:', count_hits)


if __name__ == "main":
    main(docopt.docopt(__doc__))
