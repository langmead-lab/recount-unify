#!/bin/env bash
#TODO: re-write all of this in snakemake

#params:
#1) path to top level of staging directory
#2) either samples.tsv or the datasource to use when generating sample IDs
#3) path to annotated jx collection file

set -o pipefail -o nounset -o errexit 
scripts=`perl -e '$f="'${0}'"; $f=~s/\/[^\/]+$/\//; print "$f\n";'`

#assumes this top_directory is an isolated staging area for the files (they've been moved after being found as new by the patroller)
top_directory=$1

#first create manifest files, one per loworder 2-digit grouping, to merge
#TODO: replace this with Brad's ID function rather than using samples.tsv
find ${top_directory} -name "*.jx_bed.gz" -size +0 | perl -ne 'BEGIN { open(IN,"<'${top_directory}'/samples.tsv"); %rids; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[0]; $run_acc=$f[1]; $rids{$run_acc}=$rid; } close(IN); } chomp; $f=$_; @f=split(/\//,$f); $fname=pop(@f); $run_acc=pop(@f); $loworder=pop(@f); $rid=$rids{$run_acc}; push(@{$h{$loworder}},["$f.filtered.sorted.gz",$rid,$f,$run_acc]); END { for $k (keys %h) { open(OUT,">'${top_directory}'/$k.manifest"); for $a (@{$h{$k}}) { print OUT "".join("\t",@$a)."\n"; } close(OUT); }}'

#TODO: parallelize this loop
#merge per loworder 2-digit grouping
for m in `ls ${top_directory}/*.manifest`;
do
    #filter out non-canonical splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
    #assumes ${directory} is writable
    cut -f 3 $m | perl -ne 'chomp; $f=$_; `zcat $f | egrep -e "[AG]T-A[GC]" | fgrep -v GT-AC | fgrep -v "?" | cut -f 1-7 | sort -k1,1 -k2,2n -k3,3n | gzip > $f.filtered.sorted.gz`;'
    #merge group of junctions
    python $scripts/merge.py --list-file $m --gzip | gzip > ${m}.jx_sample_files.merged.tsv.gz
    #TODO: move this to be done first (into staging) 
    cut -f 3 $m | perl -ne 'chomp; @f=split(/\//,$_); pop(@f); $f=join("\/",@f); pop(@f); $f1=join("\/",@f); $t="'${top_directory}'"; @f2=split(/\//,$t); $t=pop(@f2); $f1=~s/($t)/$1\/completed/; `mv $f $f1/`;'
done

if [ -n "${2:+1}" ]; then
	rm ${top_directory}/*.manifest
fi

#now do master merge
ls  ${top_directory}/*.jx_sample_files.merged.tsv.gz > ${top_directory}/all_jx_sample_files
#TODO need to further merge with an existing master merged file if this is an N+1 project (ongong samples being added)
#e.g. need to check if $final_dest_dir/all.jxs.merged.tsv.gz already exits and add it as the last entry to the all_jx_sample_files manifest
python $scripts/merge.py --list-file ${top_directory}/all_jx_sample_files --gzip --append-samples | gzip > ${top_directory}/all.jxs.merged.tsv.gz


#cleanup
if [ -n "${2:+1}" ]; then
	rm ${top_directory}/*.jx_sample_files.merged.tsv.gz
fi

#TODO: now annotate
