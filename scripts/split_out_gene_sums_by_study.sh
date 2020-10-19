#!/usr/bin/env bash
#assumes GNU parallel is in the path
set -o pipefail -o errexit -o nounset
#do final splits for gene sums into individual studies
s=$(dirname $0)

#e.g. sra, gtex, or tcga
comp=$1
#e.g. ERCC,SIRV,F006,G029,G026,R109" or "ERCC,SIRV,M023"
annotations=$2
#number of parallel jobs to run for splits (e.g. 40 on elephants in IDIES)
num_procs=$3

annotations_list=`echo $annotations | sed 's/,/ /g'`
for f in $annotations_list; do 
    pypy ${s}/create_gene_sums_by_study_splits.py ids.tsv ${f}.gene_sums.tsv $f gene_sums_per_study $comp 1 > gene_sums.splits.${f}.jobs
    parallel -j ${num_procs} < gene_sums.splits.${f}.jobs
done
