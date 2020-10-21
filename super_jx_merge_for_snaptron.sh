#!/usr/bin/env bash
#this will create merge, annotate, and index the complete set of unique junctions
#across multiple tranches from SRA human or mouse for Snaptron
set -exo pipefail

#requires the following to be in PATH:
#1) pypy (python 2.7 version)
#2) bgzip >=1.9
#3) tabix
#4) sqlite3

dir=$(dirname $0)

#e.g. srav3_human or srav1_mouse
#this is the prefix name of the directories where the tranches are stored,
#e.g. srav3_human5 for tranche 5 of the recount-pump recount3 output
prefix=$1
#e.g. annotated_junctions.m38.tsv.gz
annotated_jxs=$2
#compilation ID, usually 0 for SRA human, 1 for SRA mouse, 2 for GTEx, 3 for TCGA
comp_id=$3
#number of cores to use for BGZip compression
num_cpus=$4

#files like srav3_human0/all.sjs.motifs.merged.tsv
#assumes were in the directory just above the individual tranche directories
find ${prefix}? -maxdepth 1 -name "all.sjs.motifs.merged.tsv" > ${prefix}.all_jxs_motifs_merged.txt

pypy $dir/merge.py --list-file ${prefix}.all_jxs_motifs_merged.txt --append-samples --existing-sj-db "" > ${prefix}.all_jxs_merged.tsv 2> ${prefix}.all_jxs_merged.stderr

mkfifo -f import_jxs
sqlite3 ${prefix}.all_jxs_merged.sqlite < $SNAPTRON_SCHEMA_FILE

cat ${prefix}.all_jxs_merged.tsv | pypy $dir/../annotate/annotate_sjs.py --compiled-annotations $annotated_jxs --motif-correct --compilation-id $comp_id | tee import_jxs | bgzip -@ $num_cpus > ${prefix}.all_jxs_merged.annotated.bgz 2> ${prefix}.all_jxs_merged.annotated.stderr &

sqlite3 ${prefix}.all_jxs_merged.sqlite -cmd '.separator "\t"' ".import ./import_jxs intron"

#now actually index
sqlite3 ${prefix}.all_jxs_merged.sqlite < $SNAPTRON_INDEX_SCHEMA_FILE
tabix -s2 -b3 -e4 ${prefix}.all_jxs_merged.annotated.bgz

ln -fs ${prefix}.all_jxs_merged.annotated.bgz junctions.bgz
ln -fs ${prefix}.all_jxs_merged.annotated.bgz.tbi junctions.bgz.tbi
ln -fs ${prefix}.all_jxs_merged.sqlite junctions.sqlite
