#!/usr/bin/env bash
#merges generated rail_ids into original (user's) sample metadata file to create the [almost] final samples.tsv
#which is indexed by Lucene for Snaptron
set -exo pipefail

#e.g. ids.tsv
#format: <study_id>TAB<sample_id>TAB<rail_id>
sample_id_file=$1
#e.g. ids.input, or samples.input.tsv
#has to have at least 2 fields but could have many more
#format: <study_id>TAB<sample_id>....
sample_original_metadata_file=$2

#we need to join the user's original sample metadata with the newly generated rail_ids, do that using the study and sample IDs
cat $sample_original_metadata_file | perl -ne 'BEGIN { open(IN,"<'$sample_id_file'"); @ids=<IN>; close(IN); chomp(@ids); %idmap=(); map { ($study,$sample,$rid)=split(/\t/,$_); $idmap{$study."|".$sample}=$rid; } @ids; } chomp; $i++; $f=$_; @f=split(/\t/,$f); $study=shift(@f); $sample=shift(@f); $f=join("\t",@f); $f="\t$f" if(scalar @f > 0); $rail_id=$idmap{"$study|$sample"}; if(!$rail_id) { if($i==1) { print "rail_id\t$sample\t$study"."$f\n"; next; } print STDERR "failed to map $study\t$sample\t$f to a rail_id, terminating\n"; exit(-1); } print "$rail_id\t$sample\t$study"."$f\n";'
