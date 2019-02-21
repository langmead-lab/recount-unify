#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit 

#get actual jx coords; collapse ones with same coords (after sorting by coords); cut to proper input format for merge.py
cat /dev/stdin | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($o1,$o2)=split(/,/,$f[10]); $f[1]+=$o1+1; $f[2]-=$o2; print "".join("\t",@f)."\n";' | cut -f 1-3,5,6 | sort -k1,1 -k2,2n -k3,3n -k5,5 | perl -ne 'chomp; $f=$_; ($c,$s,$e,$n,$st)=split(/\t/,$f); if($pc eq $c && $ps == $s && $pe == $e) { $pn+=$n; $pst=$st; next; } else { print "".join("\t",($pc,$ps,$pe,"-1",$pn,$pst,"GT-AG"))."\n" if($pc); } ($pc,$ps,$pe,$pn,$pst)=($c,$s,$e,$n,$st); END { print "".join("\t",($pc,$ps,$pe,"-1",$pn,$pst,"GT-AG"))."\n" if($pc); }'
