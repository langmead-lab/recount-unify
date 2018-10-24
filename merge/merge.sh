#!/bin/env bash
directory=$1
#filter out non-canonical splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
#assumes ${directory} is writable
ls ${directory}/*.jx_bed | perl -ne 'chomp; $f=$_; `egrep -e "[AG]T-A[GC]"  $f | fgrep -v GT-AC | fgrep -v "?" | cut -f 1-7 | sort -k1,1 -k2,2n -k3,3n > $f.filtered.sorted`;'
#TODO: replace this with a samples metadata sourced method to create the manifest file for the merge
ls ?RR*.filtered.sorted | sort | perl -ne 'chomp; $f=$_; $f=~/_(.RP[^_]+)/; $SRP=$1; print "$f\t".$i++."\t$SRP\n";' > jx_sample_files

python merge.py sample_files > jx_sample_files.merged.tsv
