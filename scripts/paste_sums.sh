#!/bin/bash

manifest=$1
output=$2
get_ids=$3

#if not set, we're doing a final paste of all previously pasted (or copied) sample groups so no need to handle sample IDs
if [ -z ${get_ids} ]; then
    cat $manifest | perl -ne 'chomp; $f=$_; $s.=" $f"; $c++; END { if($c > 1) { print "paste $s\n"; `paste $s > '${output}'`; } else { `cp $s '${output}'`; }}'
else
#here we need to make sure we get output the correct order of samples IDs as the column header
cut -f 2,3 $manifest | perl -ne 'chomp; ($rid,$f)=split(/\t/,$_); $h.="$rid\t"; $s.=" $f"; $c++; END { $h=~s/\t$//; open(OUT,">'${output}'"); print OUT "$h\n"; close(OUT); if($c > 1) { print "paste $s\n"; `paste $s >> '${output}'`; } else { `cat $s >> '${output}'`; }}'
fi
