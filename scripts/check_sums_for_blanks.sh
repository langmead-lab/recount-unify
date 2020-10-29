#!/usr/bin/bash
#checks for empty strings where numbers (>=0) should be
#in the base disjoint and rejoined exon sums files

dir=$(dirname $0)
#run check_for_blanks
tranche_dir=$1

cd $tranche_dir
num_fields=`pigz --stdout -p2 -d all.exon_bw_count.pasted.gz | head -2 | tail -n1 | tr \\\\t \\\\n | wc -l`
/usr/bin/time -v pcat all.exon_bw_count.pasted.gz 2> check_for_blanks.timing2 2>&1 | ${dir}/check_for_blanks -n $num_fields > check_for_blanks.blanks2 2> check_for_blanks.run2

num_fields=`pigz --stdout -p2 -d all.exon_counts.rejoined.tsv.gz | head -2 | tail -n1 | tr \\\\t \\\\n | wc -l`
/usr/bin/time -v pcat all.exon_counts.rejoined.tsv.gz 2> check_for_blanks.exons_rejoined.timing2 2>&1 | ${dir}/check_for_blanks -n $num_fields > check_for_blanks.exons_rejoined.blanks 2> check_for_blanks.exons_rejoined.run
