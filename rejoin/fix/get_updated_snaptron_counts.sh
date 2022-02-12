#!/usr/bin/env bash
set -exo pipefail
annot="G026"

study=$1
pushd $study

updated_gene_sums=sra.gene_sums.${study}.${annot}.gz
set +o pipefail
pcat $updated_gene_sums | tail -n+3 | head -1 | cut -f 2- | tr '\t' '\n' > ${updated_gene_sums}.samples
set -o pipefail

/usr/bin/time -v pcat $updated_gene_sums | tail -n+4 | fgrep -f /datascope/recount03/rejoin_fix/disjoint2exons2genes.fix.sorted.bed.${annot}_genes.tabs  | python3 /datascope/recount03/recount-unify/rejoin/fix/find_updated_gene_counts.py ${updated_gene_sums}.samples <(fgrep $'\t'"$study"$'\t' ../../recount3.sra.samples.tsv.cut | cut -f 1-3) | LC_ALL=C sort -k1,1 | cut -f2- > ${annot}.snaptron

echo "DONE"
popd
