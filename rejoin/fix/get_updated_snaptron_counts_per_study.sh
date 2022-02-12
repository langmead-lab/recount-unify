#!/usr/bin/env bash
set -exo pipefail
dir=$(dirname $0)
ddir="/datascope/recount03/rejoin_fix"
annot="G026"
genes="$ddir/disjoint2exons2genes.fix.sorted.bed.${annot}_genes.tabs"
ids="$ddir/recount3.sra.samples.tsv.cut"

study=$1

cwd=$(basename $PWD)
updated_gene_sums=sra.gene_sums.${study}.${annot}.gz
if [[ $cwd == "tcga." ]]; then
    ids="$ddir/tcgav2.samples.tsv.cut"
    updated_gene_sums=tcga.gene_sums.${study}.${annot}.gz
fi
if [[ $cwd == "gtex." ]]; then
    ids="$ddir/gtexv2.samples.tsv.cut"
    updated_gene_sums=gtex.gene_sums.${study}.${annot}.gz
fi
    
pushd $study

set +o pipefail
pcat $updated_gene_sums | tail -n+3 | head -1 | cut -f 2- | tr '\t' '\n' > ${updated_gene_sums}.samples
set -o pipefail

if [[ $cwd == "gtex." ]]; then
    /usr/bin/time -v pcat $updated_gene_sums | tail -n+4 | fgrep -f $genes | python3 $dir/find_updated_gene_counts.py ${updated_gene_sums}.samples <(fgrep $'\t'"$study"$'\t' $ids | cut -f 1,3-4) | LC_ALL=C sort -k1,1 | cut -f2- > ${annot}.snaptron
else
    /usr/bin/time -v pcat $updated_gene_sums | tail -n+4 | fgrep -f $genes | python3 $dir/find_updated_gene_counts.py ${updated_gene_sums}.samples <(fgrep $'\t'"$study"$'\t' $ids | cut -f 1-3) | LC_ALL=C sort -k1,1 | cut -f2- > ${annot}.snaptron
fi

echo "DONE"
popd
