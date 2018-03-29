#!/bin/bash
indir=$1

for f in $(ls $indir); do locus=` echo $f | sed 's/.fna//' | sed 's/.haps//'`; grep '>' ${indir}$f | sed 's/>//' | sed 's/_.*//' >> ${locus}_ids_unsorted.txt; done
for f in *unsorted.txt; do new=`echo $f | sed 's/_unsorted//'`; tail -n +2 $f | sort | uniq > $new; done
rm *_unsorted.txt
for f in *.txt; do wc -l $f | sed 's/_ids.txt//' >> count_individuals_per_locus.txt; done
