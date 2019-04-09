#!/bin/bash
set -o pipefail -o nounset -o errexit 

date=`date +%s`

search_dir=$1
destination_dir=$2

#make top level (study) pre-staging directory
mkdir -p ${destination_dir}

find -L $search_dir -name "*.done" | sort | perl -ne 'BEGIN { $dest_dir="'${destination_dir}'"; } chomp; $fp=$_; $fp=~s/\.done$//; @f=split(/\//,$fp); $f=pop(@f); $acc=pop(@f); $acc_loworder=pop(@f); $study=pop(@f); $loworder=pop(@f); ($proj,$input_id,$attempt)=split(/_/,$f); $f_wo_attempt_num=$f; $f_wo_attempt_num=~s/(\d+)$//; $attempt_num=$1; print "$study\t$f\t$loworder\t$f_wo_attempt_num\t$attempt_num\t$acc_loworder\n"; if(!$h{$study.$f_wo_attempt_num} || $h{$study.$f_wo_attempt_num}->[0] > $attempt_num) { $h{$study.$f_wo_attempt_num}=[$attempt_num,$fp,$loworder,$study,$acc_loworder,$acc]; } END { for $v (values %h) { ($attempt_num,$fp,$loworder,$study,$acc_loworder,$acc)=@$v; `mkdir -p $dest_dir/$loworder/$study/$acc_loworder/$acc`; print "$fp\n"; `ln -s $fp $dest_dir/$loworder/$study/$acc_loworder/$acc/`; }}' > ${destination_dir}/attempts.moved.${date}
