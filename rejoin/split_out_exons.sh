#!/usr/bin/env bash
set -o pipefail -o errexit -o nounset
#use this to do the following to get per-study, per-annotation exon sums from rejoined exons sums in a tranche/project
#1) Cuts out all the columns which are samples in the study from the rejoined exons file
#2) Further splits the cut file into one or more annotation sums files
#3) Compresses the resulting per-annotation files

#2020-04-29 15:58:42
date=`date "+%Y-%m-%d %H:%M:%S"`

#src='sra'
src='gtex'

dir=$(dirname $0)

#e.g. G026,G029,R109,ERCC,SIRV
annotations=$1
#e.g. all.exon_counts.rejoined.tsv.gz.coords.bitmasks2 or
#exon_counts.bitmasks.tsv
#produced by map_exon_sum_rows_to_annotations.py
#e.g. exon_counts.bitmasks.tsv
annotation_bitmasks=$2
#e.g. 1709834 for human recount3
num_exons=$3
#e.g. ERP001942
study=$4
study_file=$5
target_dir=$6
#e.g. exon_counts.bitmasks.coords
coords_file=$7
src=$8

mkdir -p $target_dir

#split out annotations
#assumes there's a header in each per-study file
pcat $study_file | cut -f 2- | $dir/split_exons -a $annotations -b $annotation_bitmasks -c $coords_file -n $num_exons -p $study -h

#name format e.g. sra.exon_sums.SRP186687.G026.gz
annotations_list=`echo $annotations | sed 's/,/ /g'`
for a in $annotations_list; do cat <(echo "##annotation=$a") <(echo "##date.generated=$date") <(paste <(echo -n "chromosome|start_1base|end_1base|strand") <(pcat $study_file | head -1 | cut -f 2-)) <(cat ${a}.${study}.tsv) | pigz --fast -p1 > $target_dir/${src}.exon_sums.${study}.${a}.gz ; rm ${a}.${study}.tsv ; done
