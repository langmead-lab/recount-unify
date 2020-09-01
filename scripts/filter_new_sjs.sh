#!/bin/bash
set -o pipefail -o errexit 

manifest=$1

#[optional] number of threads to use when pigz'ing output, default 1
threads=$2
if [[ -z $threads ]]; then
	threads=1
fi

#reformat splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
#assumes $f's directory is writable

#ensure sort order is the same across all
export LC_ALL=C
#motifs are filled in later rather than here
for f in `cut -f 3 ${manifest}`; do
    if [[ $threads -gt 0 ]]; then
        zstd -cd $f 2>/dev/null | perl -ne 'chomp; $f=$_; ($chr,$start,$end,$strand,$motif,$annotated,$num_uniques,$num_multis,$max_overhang)=split(/\t/,$f); $total_cov=$num_uniques+$num_multis; $strand=~tr/012/\?\+\-/; print "$chr\t$start\t$end\t$num_uniques\t$total_cov\t$strand\tCC-CC\n";' | sort -k1,1 -k2,2n -k3,3n -u | pigz --fast -p $threads > ${f}.filtered.sorted.gz;
    else
        zstd -cd $f 2>/dev/null | perl -ne 'chomp; $f=$_; ($chr,$start,$end,$strand,$motif,$annotated,$num_uniques,$num_multis,$max_overhang)=split(/\t/,$f); $total_cov=$num_uniques+$num_multis; $strand=~tr/012/\?\+\-/; print "$chr\t$start\t$end\t$num_uniques\t$total_cov\t$strand\tCC-CC\n";' | sort -k1,1 -k2,2n -k3,3n -u > ${f}.filtered.sorted
    fi
done
touch ${manifest}.filtered
