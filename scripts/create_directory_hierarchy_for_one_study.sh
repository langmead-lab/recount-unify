#!/bin/bash
set -o pipefail -o nounset -o errexit 

#e.g. ids.tsv
ids_file=$1

#expects full path, otherwise the symlinking won't work
input_dir=$2

#this can be a relative path
#this is not the final links directory, rather a temporary one which will further be linked from
link_dir=$3

cat <(cut -f 1 $ids_file | head -1) <(ls $input_dir) | perl -ne 'BEGIN { $input_dir="'$input_dir'"; $link_dir="'$link_dir'"; } chomp; if(!$study) { $study=$_; next; } $sample_dir=$_; $study=~/(..)$/; $lo1=$1; $sample_dir=~/^(.*(..))_att\d+$/; $sample=$1; $lo2=$2; `mkdir -p $link_dir/$lo1/$study/$lo2/$sample`; `ln -fs $input_dir/$sample_dir $link_dir/$lo1/$study/$lo2/$sample/$sample_dir`; `touch $link_dir/$lo1/$study/$lo2/$sample/$sample_dir.done`;'
