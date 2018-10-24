#!/bin/env bash
directory=$1
#filter out non-canonical splice junctions and sort by coordinates, while cutting out extraneous fields, and ? stranded junctions
#assumes ${directory} is writable
ls ${directory}/*.jx_bed | perl -ne 'chomp; $f=$_; `egrep -e "[AG]T-A[GC]"  $f | fgrep -v GT-AC | fgrep -v "?" | cut -f 1-7 | sort -k1,1 -k2,2n -k3,3n > $f.filtered.sorted`;'
#TODO: replace this with a samples metadata sourced method to create the manifest file for the merge
ls ?RR*.filtered.sorted | sort | perl -ne 'chomp; $f=$_; $f=~/_(.RP[^_]+)/; $SRP=$1; print "$f\t".$i++."\t$SRP\n";' > jx_sample_files

python merge.py sample_files > jx_sample_files.merged.tsv
....multiple merge steps later.....
#after all studies have been merged, generate per-study jx files (both coverage and BED) for recount
zcat jx_sample_files.merged.tsv.gz | perl -ne 'BEGIN { open(IN,"<samples.tsv"); %h; %fhs; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[0]; $study=$f[4]; $h{$rid}=$study; if(!$fhs{$study}) { $fhs{$study}->[0]=IO::File->new("|gzip >$study.junction_coverage.tsv.gz"); $fhs{$study}->[1]=IO::File->new("|gzip >$study.junction_id_with_transcripts.bed.gz"); }} close(IN); $jx_id=-1; } chomp; $jx_id++; ($chrm,$start,$end,$strand,$samples_covs)=split(/\t/,$_); $samples_covs=~s/^,//; %studies_seen=(); for $sample (split(/,/,$samples_covs)) { ($sid,$c)=split(/:/,$sample); $study=$h{$sid}; ($covF,$bedF)=@{$fhs{$study}}; print $covF "$jx_id\t$sid\t$c\n"; if(!$studies_seen{$study}) { print $bedF "$chrm\t$start\t$end\t$jx_id|D:NA|A:NA|J:NA\t1000\t$strand\n"; $studies_seen{$study}=1; } }'
