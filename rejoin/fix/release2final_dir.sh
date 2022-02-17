#!/usr/bin/env bash
set -exo pipefail
dir=$(dirname $0)

study_listF=$1
src=$2

ddir=/datascope/recount03/release/human/data_sources/$src/gene_sums_postrejoinfix

#for study_listF in `ls sra.studies??`; do
sdir="DONES/"$(echo "$study_listF" | sed 's#studies##')
pushd $sdir
for s in `cat ../../$study_listF`; do
    lo=${$s: -2}
    target_dir=$ddir/$lo/$s
    mkdir -p $target_dir
    for annot in G026 G029 R109 F006; do
        ln -f $s/${src}.gene_sums.${s}.${annot}.gz $target_dir/${src}.gene_sums.${s}.${annot}.gz
    done
done
