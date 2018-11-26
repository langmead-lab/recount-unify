#!/bin/bash
set -o pipefail -o errexit 

manifest=$1
output=$2
dont_get_ids=$3
existing_sums=$4

if [ -n "$dont_get_ids" ]; then
#if set, we're doing a final paste of all previously pasted (or copied) sample groups so no need to handle sample IDs
    cat $manifest | perl -ne 'chomp; $f=$_; $s.=" $f"; $c++; END { if($c > 1) { print "paste $s\n"; `paste $s | gzip > '${output}'`; } else { `cp $f '${output}'`; }}'
else
#here we need to make sure we output the correct order of samples IDs as the column header
    cut -f 2,3 $manifest | perl -ne 'chomp; ($rid,$f)=split(/\t/,$_); $h.="$rid\t"; $s.=" $f.unc"; $c++; END { $h=~s/\t$//; open(OUT,">'${output}'"); print OUT "$h\n"; close(OUT); if($c > 1) { print "paste $s\n"; `paste $s >> '${output}'`; } else { `cat $f.unc >> '${output}'`; }}'
fi

if [ -n "$existing_sums" ]; then
    mv ${output} ${output}.pre_existing
    paste <(zcat ${existing_sums}) <(zcat ${output}.pre_existing) | gzip > ${output}
fi
