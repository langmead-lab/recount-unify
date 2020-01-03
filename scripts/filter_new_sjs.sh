#!/bin/bash
set -o pipefail -o nounset -o errexit 

#reformat splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
#assumes $f's directory is writable

#motifs are filled in later rather than here
for f in `cut -f 3 ${1}`; do
    zstd -cd $f 2>/dev/null | perl -ne 'chomp; $f=$_; ($chr,$start,$end,$strand,$motif,$annotated,$num_uniques,$num_multis,$max_overhang)=split(/\t/,$f); $total_cov=$num_uniques+$num_multis; $strand=~tr/012/\?\+\-/; print "$chr\t$start\t$end\t$num_multis\t$total_cov\t$strand\tCC-CC\n";' | sort -k1,1 -k2,2n -k3,3n -u | gzip > ${f}.filtered.sorted.gz;
done
touch ${1}.filtered
