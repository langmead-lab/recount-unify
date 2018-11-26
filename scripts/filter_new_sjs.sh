#!/bin/bash
set -o pipefail -o nounset -o errexit 

#filter out non-canonical splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
#assumes $f's directory is writable
for f in `cut -f 3 ${1}`; do
    zstd -cd $f 2>/dev/null | perl -ne 'BEGIN { %MOTIF_MAP=(1=>"GT-AG",2=>"GT-AG",3=>"GC-AG",4=>"GC-AG",5=>"AT-AC",6=>"AT-AC"); } chomp; $f=$_; ($chr,$start,$end,$strand,$motif,$annotated,$num_uniques,$num_multis,$max_overhand)=split(/\t/,$f); next if(!$MOTIF_MAP{$motif} || $strand == 0 || ($motif=~/[246]/ && $strand == 1) || ($motif=~/[135]/ && $strand == 2)); $total_cov=$num_uniques+$num_multis; $strand=~tr/12/\+\-/; print "$chr\t$start\t$end\t$num_multis\t$total_cov\t$strand\t".$MOTIF_MAP{$motif}."\n";' | sort -k1,1 -k2,2n -k3,3n | gzip > ${f}.filtered.sorted.gz;
done
touch ${1}.filtered
