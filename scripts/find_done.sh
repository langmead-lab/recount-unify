#!/bin/bash
set -o pipefail -o nounset -o errexit 

date=`date +%s`

script_dir=$(dirname $0)

search_dir=$1
destination_dir=$2


#make top level (study) pre-staging directory
mkdir -p ${destination_dir}
find -L $search_dir -name "sra_human_v3*" -type d | fgrep "att" > ${destination_dir}/all_attempt_dirs

#TODO: minor bug here while this will favor eariler attempts over later ones, except where the earliest extant attempt is 0 *and* it comes first in the list, 
#then it will favor the next earliest
#keeping it for now to keep all sub-compilations consistent
find -L $search_dir -name "*.done" | sort | perl -ne 'BEGIN { $dest_dir="'${destination_dir}'"; } chomp; $fp=$_; $fp=~s/\.done$//; @f=split(/\//,$fp); $f=pop(@f); $fp2=join("/",@f); $acc=pop(@f); $acc_loworder=pop(@f); $study=pop(@f); $loworder=pop(@f); ($proj,$input_id,$attempt)=split(/_/,$f); $f_wo_attempt_num=$f; $f_wo_attempt_num=~s/(\d+)$//; $attempt_num=$1; print "$study\t$f\t$loworder\t$f_wo_attempt_num\t$attempt_num\t$acc_loworder\n"; $fp=$fp2; if(!$h{$study.$f_wo_attempt_num} || $h{$study.$f_wo_attempt_num}->[0] > $attempt_num) { $h{$study.$f_wo_attempt_num}=[$attempt_num,$fp,$loworder,$study,$acc_loworder,$acc]; } END { for $v (values %h) { ($attempt_num,$fp,$loworder,$study,$acc_loworder,$acc)=@$v; `mkdir -p $dest_dir/$loworder/$study/$acc_loworder`; print "$fp\n"; `ln -fs $fp $dest_dir/$loworder/$study/$acc_loworder/`; }}' > ${destination_dir}/attempts.moved.${date}

mv ${destination_dir}/attempts.moved.${date} ${destination_dir}/attempts.moved.complete

#now remove redundant attempts
fgrep "/" ${destination_dir}/attempts.moved.complete > ${destination_dir}/completed_attempt_dirs
python ${script_dir}/filter_one_file_by_another.py -f ${destination_dir}/completed_attempt_dirs -t ${destination_dir}/all_attempt_dirs -p -w -c 0 -n > ${destination_dir}/all_attempt_dirs.to_delete
cat ${destination_dir}/all_attempt_dirs.to_delete | perl -ne 'chomp; `rm -rf $_`;'
