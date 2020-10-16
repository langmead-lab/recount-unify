#!/usr/bin/env bash
#assumes GNU parallel is in the path
set -o pipefail -o errexit -o nounset
#do final splits for exon sums into individual studies
s=$(dirname $0)

#e.g. sra, gtex, or tcga
src=$1
#e.g. G026,G029,R109,F006,ERCC,SIRV or M023,ERCC,SIRV
annotations=$2
#e.g. 1709834 for human recount3 run
num_exons=$3
#e.g. exon_counts.bitmasks.tsv, contains lines like:
#"000000000000001000      ERCC-00002      ERCC-00002      0       1061    1061    +"
bitmasks_file=$4
#e.g. exon_counts.bitmasks.coords, contains likes like: 
#"ERCC-00002|1|1061|+"
bitmasks_coords_file=$5
#e.g. all.exon_counts.rejoined.tsv.gz
rejoin_exon_sums=$6
#number of parallel jobs to run for splits (e.g. 40 on elephants in IDIES)
num_procs=$7

exon_split_by_study_dir='exons_split_by_study_temp'
mkdir -p $exon_split_by_study_dir

#first do split-by-study (not annotation yet)
pypy $s/create_exon_sums_by_study_splits.py ids.tsv $rejoin_exon_sums ${rejoin_exon_sums}.accession_header $exon_split_by_study_dir > exon_sums.study_splits.jobs

parallel -j $num_procs < exon_sums.study_splits.jobs

#now do final splits into study.annotation files
mkdir -p exon_annotation_split_runs
ls $exon_split_by_study_dir/*.gz | perl -ne 'chomp; $f=$_; ($dir,$study)=split(/\//,$f); $study=~s/\.gz//; @f2=split(//,$study); $lo=pop(@f2); $lo=pop(@f2).$lo; $dir="exon_sums_per_study/$lo/$study"; `mkdir -p $dir`; print "/bin/bash -x '$s'/split_out_exons.sh \"'$annotations'\" '$bitmasks_file' '$num_exons' $study '$exon_split_by_study_dir'/$study.gz $dir '$bitmasks_coords_file' '$src' > exon_annotation_split_runs/$study.run 2>&1\n";' > exon_sums.annotation_splits.jobs

parallel -j $num_procs < exon_sums.annotation_splits.jobs
