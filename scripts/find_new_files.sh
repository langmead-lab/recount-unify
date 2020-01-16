#!/bin/bash
set -o pipefail -o nounset -o errexit 


DELIM='!';

#top level of incoming dir
#(where fully done analyses were copied/symlinked from original set of attempts)
search_dir=$1
#map between Rail generated IDs and externall accessions/UUIDs, format: rail_id<TAB>accession
sample_ids_file=$2
destination_dir=$3
#e.g. "unique.exon_bw_count" or "sj" or "all.exon_bw_count"
analysis_type=$4
#.e.g ".zst"
additional_suffix=$5
#optional flag, if set, we're doing *per-study* not per-study-loworder digits
per_study=$6

#now uses both study loworder digits *and* run loworder digits in manifest file names
if [[ -z $per_study ]]; then
	find -L $search_dir -name "*.${analysis_type}${additional_suffix}" -size +0 | perl -ne 'BEGIN { open(IN,"<'${sample_ids_file}'"); %rids; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[2]; $run_acc=$f[1]; $rids{$run_acc}=$rid; } close(IN); } chomp; $f=$_; @f=split(/\//,$f); $fname=pop(@f); `pushd '${destination_dir}' && ln -fs ../$f && popd`; $fname=~/^([^!]+)!/; $run_acc=$1; $attempt=pop(@f); $acc=pop(@f); $acc_loworder=pop(@f); $study=pop(@f); $study_loworder=pop(@f); $rid=$rids{$run_acc}; push(@{$h{$study_loworder.".".$acc_loworder}},["'${destination_dir}'/$fname.filtered.sorted.gz",$rid,"'${destination_dir}'/$fname",$run_acc]); END { open(ALL_OUT,">'${destination_dir}'/'${analysis_type}'.groups.manifest"); for $k (keys %h) { open(OUT,">'${destination_dir}'/'${analysis_type}'.$k.manifest"); for $a (@{$h{$k}}) { print ALL_OUT "".$a->[3]."\t'${analysis_type}'.$k.manifest\n"; print OUT "".join("\t",@$a)."\n"; } close(OUT); } close(ALL_OUT); }'

else
	find -L $search_dir -name "*.${analysis_type}${additional_suffix}" -size +0 | perl -ne 'BEGIN { open(IN,"<'${sample_ids_file}'"); %rids; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[2]; $run_acc=$f[1]; $rids{$run_acc}=$rid; } close(IN); } chomp; $f=$_; @f=split(/\//,$f); $fname=pop(@f); `pushd '${destination_dir}' && ln -fs ../$f && popd`; $fname=~/^([^!]+)!/; $run_acc=$1; $attempt=pop(@f); $acc=pop(@f); $acc_loworder=pop(@f); $study=pop(@f); $rid=$rids{$run_acc}; push(@{$h{$study.".".$acc_loworder}},["'${destination_dir}'/$fname.filtered.sorted.gz",$rid,"'${destination_dir}'/$fname",$run_acc]); END { open(ALL_OUT,">'${destination_dir}'/'${analysis_type}'.groups.manifest"); for $k (keys %h) { open(OUT,">'${destination_dir}'/'${analysis_type}'.$k.manifest"); for $a (@{$h{$k}}) { print ALL_OUT "".$a->[3]."\t'${analysis_type}'.$k.manifest\n"; print OUT "".join("\t",@$a)."\n"; } close(OUT); } close(ALL_OUT); }'
fi
