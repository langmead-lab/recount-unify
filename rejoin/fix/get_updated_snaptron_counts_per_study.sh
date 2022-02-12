#!/usr/bin/env bash
set -exo pipefail
dir=$(dirname $0)
ddir="/datascope/recount03/rejoin_fix"
annot="G026"
genes="$ddir/disjoint2exons2genes.fix.sorted.bed.${annot}_genes.tabs"
ids="$ddir/recount3.sra.samples.tsv.cut"

study=$1
pushd $study

updated_gene_sums=sra.gene_sums.${study}.${annot}.gz
set +o pipefail
pcat $updated_gene_sums | tail -n+3 | head -1 | cut -f 2- | tr '\t' '\n' > ${updated_gene_sums}.samples
set -o pipefail

/usr/bin/time -v pcat $updated_gene_sums | tail -n+4 | fgrep -f $genes | python3 $dir/find_updated_gene_counts.py ${updated_gene_sums}.samples <(fgrep $'\t'"$study"$'\t' $ids | cut -f 1-3) | LC_ALL=C sort -k1,1 | cut -f2- > ${annot}.snaptron

echo "DONE"
popd
