#!/bin/bash
set -o pipefail -o nounset -o errexit 

search_dir=$1
#map between Rail generated IDs and externall accessions/UUIDs, format: rail_id<TAB>accession
sample_ids_file=$2
destination_dir=$3
analysis_type=$4
additional_suffix=$5

find $search_dir -name "*.${analysis_type}${additional_suffix}" -size +0 | perl -ne 'BEGIN { open(IN,"<'${sample_ids_file}'"); %rids; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[0]; $run_acc=$f[1]; $rids{$run_acc}=$rid; } close(IN); } chomp; $f=$_; `rsync -av $f '${destination_dir}'/`; @f=split(/\//,$f); $fname=pop(@f); $fname=~/^([^_]+)_/; $run_acc=$1; $attempt=pop(@f); $loworder=pop(@f); $rid=$rids{$run_acc}; push(@{$h{$loworder}},["'${destination_dir}'/$fname.filtered.sorted.gz",$rid,"'${destination_dir}'/$fname",$run_acc]); END { open(ALL_OUT,">'${destination_dir}'/'${analysis_type}'.groups.manifest"); for $k (keys %h) { open(OUT,">'${destination_dir}'/'${analysis_type}'.$k.manifest"); for $a (@{$h{$k}}) { print ALL_OUT "".$a->[3]."\t'${analysis_type}'.$k.manifest\n"; print OUT "".join("\t",@$a)."\n"; } close(OUT); } close(ALL_OUT); }'
