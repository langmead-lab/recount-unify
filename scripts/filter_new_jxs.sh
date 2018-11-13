#!/bin/bash
#filter out non-canonical splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
#assumes ${directory} is writable
cut -f 3 $1 | perl -ne 'chomp; $f=$_; `zcat $f | egrep -e "[AG]T-A[GC]" | fgrep -v GT-AC | fgrep -v "?" | cut -f 1-7 | sort -k1,1 -k2,2n -k3,3n | gzip > $f.filtered.sorted.gz`;'
touch ${1}.filtered
