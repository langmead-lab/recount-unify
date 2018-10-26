#!/bin/env bash
set -o pipefail -o nounset -o errexit 
scripts=`perl -e '$f="'${0}'"; $f=~s/\/[^\/]+$/\//; print "$f\n";'`

top_directory=$1

#first create manifest files, one per loworder 2-digit grouping, to merge
#ls ${top_directory}/ | fgrep manifest | perl -ne 'chomp; `rm $_`;'
find ${top_directory} -name "*.jx_bed.gz" -size +0 | perl -ne 'BEGIN { open(IN,"<'${top_directory}'/samples.tsv"); %rids; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[0]; $run_acc=$f[1]; $rids{$run_acc}=$rid; } close(IN); } chomp; $f=$_; @f=split(/\//,$f); $fname=pop(@f); $run_acc=pop(@f); $loworder=pop(@f); $rid=$rids{$run_acc}; push(@{$h{$loworder}},["$f.filtered.sorted.gz",$rid,$f]); END { for $k (keys %h) { open(OUT,">'${top_directory}'/$k.manifest"); for $a (@{$h{$k}}) { print OUT "".join("\t",@$a)."\n"; } close(OUT); }}'

#ls ${top_directory}/ | fgrep jx_sample_files.merged.tsv.gz | perl -ne 'chomp; `rm $_`;'
#TODO: parallelize this loop
#merge per loworder 2-digit grouping
for m in `ls ${top_directory}/*.manifest`;
do
    #filter out non-canonical splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
    #assumes ${directory} is writable
    cut -f 3 $m | perl -ne 'chomp; $f=$_; `zcat $f | egrep -e "[AG]T-A[GC]" | fgrep -v GT-AC | fgrep -v "?" | cut -f 1-7 | sort -k1,1 -k2,2n -k3,3n | gzip > $f.filtered.sorted.gz`;'
    #TODO: replace this with a samples metadata sourced method to create the manifest file for the merge
    #cut -f 1 $m | sort | perl -ne 'chomp; $f=$_.".filtered.sorted"; $f=~/_(.RP[^_]+)/; $SRP=$1; print "$f\t".$i++."\t$SRP\n";' > jx_sample_files
    python $scripts/merge.py --list-file $m --gzip | gzip > ${m}.jx_sample_files.merged.tsv.gz
done

if [ -n "${2:+1}" ]; then
	rm ${top_directory}/*.manifest
fi

#now do master merge
ls  ${top_directory}/*.jx_sample_files.merged.tsv.gz > ${top_directory}/all_jx_sample_files
#TODO: change merge.py to handle samples at this level differently
python $scripts/merge.py --list-file ${top_directory}/all_jx_sample_files --gzip --append-samples | gzip > ${top_directory}/all.jxs.merged.tsv.gz
 
#after all studies have been merged, generate per-study jx files (both coverage and BED) for recount
zcat ${top_directory}/all.jxs.merged.tsv.gz | perl -ne 'BEGIN { $suffix1="junction_id_with_transcripts.bed.gz"; $suffix2="junction_coverage.tsv.gz"; open(IN,"<'${top_directory}'/samples.tsv"); %h; %fhs; %used; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[0]; $study=$f[4]; $h{$rid}=$study; if(!$fhs{$study}) { $fhs{$study}->[0]=IO::File->new("|gzip >'${top_directory}'/$study.$suffix1"); $fhs{$study}->[1]=IO::File->new("|gzip >'${top_directory}'/$study.$suffix2"); }} close(IN); $jx_id=-1; } chomp; $jx_id++; ($chrm,$start,$end,$strand,$samples_covs)=split(/\t/,$_); $samples_covs=~s/^,//; %studies_seen=(); for $sample (split(/,/,$samples_covs)) { ($sid,$c)=split(/:/,$sample); $study=$h{$sid}; $used{$study}=1; ($covF,$bedF)=@{$fhs{$study}}; print $covF "$jx_id\t$sid\t$c\n"; if(!$studies_seen{$study}) { print $bedF "$chrm\t$start\t$end\t$jx_id|D:NA|A:NA|J:NA\t1000\t$strand\n"; $studies_seen{$study}=1; }} END { for $study (values %h) { if(!$used{$study}) { $f1="'${top_directory}'/$study.$suffix1"; $f2="'${top_directory}'/$study.$suffix2"; `rm $f1` if(-e $f1); `rm $f2` if(-e $f2); $used{$study}=1;}}}'

if [ -n "${2:+1}" ]; then
	rm ${top_directory}/*.jx_sample_files.merged.tsv.gz
fi
