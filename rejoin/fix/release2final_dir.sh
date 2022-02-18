#!/usr/bin/env bash
set -eo pipefail
dir=$(dirname $0)

study_listF=$1
src=$2

originald=/datascope/recount03/release/human/data_sources/$src/gene_sums
ddir=/datascope/recount03/release/human/data_sources/$src/gene_sums_postrejoinfix

#for study_listF in `ls sra.studies??`; do
sdir="DONES/"$(echo "$study_listF" | sed 's#studies##')
pushd $sdir
for s in `cat ../../$study_listF`; do
    lo=${s: -2}
    target_dir=$ddir/$lo/$s
    #echo "mkdir -p $target_dir"
    mkdir -p $target_dir
    for annot in G026 G029 R109 F006; do
        #echo "ln -f $s/${src}.gene_sums.${s}.${annot}.gz $target_dir/${src}.gene_sums.${s}.${annot}.gz"
        ln -f $s/${src}.gene_sums.${s}.${annot}.gz $target_dir/${src}.gene_sums.${s}.${annot}.gz
    done
    original_target_dir=$originald/$lo/$s
    ln -f $original_target_dir/${src}.gene_sums.${s}.ERCC.gz $target_dir/${src}.gene_sums.${s}.ERCC.gz
    ln -f $original_target_dir/${src}.gene_sums.${s}.SIRV.gz $target_dir/${src}.gene_sums.${s}.SIRV.gz
done
popd
