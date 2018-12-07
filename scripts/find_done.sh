#!/bin/bash
set -o pipefail -o nounset -o errexit 

#date=`date +%Y%m%d.%s`
date=`date +%s`

search_dir=$1
destination_dir=$2
#samples=$3

#make top level (study) pre-staging directory
perl -e 'if(! -d "'${destination_dir}'") { `mkdir -p '${destination_dir}'`; }'

#make loworder (per-study) pre-staging directories
#UPDATE: dont do this, leads to errors downstream when empty directories are picked up
#perl -e 'if(! -d "'${destination_dir}'/99") { for $i (0..99) { $i1=$i; $i1="0".$i1 if(length($i) == 1); `mkdir '${destination_dir}'/$i1`;} }'

#now 1) rsync the latest attempt directory which has a .done file and 2) delete the rest of the attemps/remaining
#find $search_dir -name "*.done" | cut -d'.' -f 1 | sort | perl -ne 'BEGIN { $dest_dir="'${destination_dir}'"; } chomp; $fp=$_; @f=split(/\//,$fp); $f=pop(@f); $loworder=pop(@f); ($proj,$input_id,$attempt)=split(/_/,$f); $f_wo_attempt_num=$f; $f_wo_attempt_num=~s/\d+$//; $h{$f_wo_attempt_num}=[$fp,$loworder]; END { for $v (values %h) { ($fp,$loworder)=@$v; if(! -d "$dest_dir/$loworder") { `mkdir -p $dest_dir/$loworder`; } print "$fp\n"; `rsync -av $fp $dest_dir/$loworder/`; $fp_no_attempt=$fp; $fp_no_attempt=~s/\d+$//; `rm -rf $fp_no_attempt*`; }}' > ${destination_dir}/attempts.moved.${date}
#dont do the rm yet
#also do hard links instead of rsyncs
#find -L $search_dir -name "*.done" | sort | perl -ne 'BEGIN { open(IN,"<"'${samples}'"); while($line=<IN>) { chomp($line); ($rid,$srr,$srp)=split(/\t/,$line); $h{$srr}=$srp; } close(IN); $dest_dir="'${destination_dir}'"; } chomp; $fp=$_; $fp=~s/\.done$//; @f=split(/\//,$fp); $f=pop(@f); $study=pop(@f); $loworder=pop(@f); ($proj,$input_id,$attempt)=split(/_/,$f); $f_wo_attempt_num=$f; $f_wo_attempt_num=~s/(\d+)$//; $attempt_num=$1; print "$study\t$f\t$loworder\t$f_wo_attempt_num\t$attempt_num\n"; if(!$h{$study.$f_wo_attempt_num} || $h{$study.$f_wo_attempt_num}->[0] > $attempt_num) { $h{$study.$f_wo_attempt_num}=[$attempt_num,$fp,$loworder,$study]; } END { for $v (values %h) { ($attempt_num,$fp,$loworder,$study)=@$v; if(! -d "$dest_dir/$loworder/$study") { `mkdir -p $dest_dir/$loworder/$study`; } print "$fp\n"; `ln $fp $dest_dir/$loworder/$study/`; }}' > ${destination_dir}/attempts.moved.${date}
find -L $search_dir -name "*.done" | sort | perl -ne 'BEGIN { $dest_dir="'${destination_dir}'"; } chomp; $fp=$_; $fp=~s/\.done$//; @f=split(/\//,$fp); $f=pop(@f); $study=pop(@f); $loworder=pop(@f); ($proj,$input_id,$attempt)=split(/_/,$f); $f_wo_attempt_num=$f; $f_wo_attempt_num=~s/(\d+)$//; $attempt_num=$1; print "$study\t$f\t$loworder\t$f_wo_attempt_num\t$attempt_num\n"; if(!$h{$study.$f_wo_attempt_num} || $h{$study.$f_wo_attempt_num}->[0] > $attempt_num) { $h{$study.$f_wo_attempt_num}=[$attempt_num,$fp,$loworder,$study]; } END { for $v (values %h) { ($attempt_num,$fp,$loworder,$study)=@$v; if(! -d "$dest_dir/$loworder/$study") { `mkdir -p $dest_dir/$loworder/$study`; } print "$fp\n"; `ln $fp $dest_dir/$loworder/$study/`; }}' > ${destination_dir}/attempts.moved.${date}
