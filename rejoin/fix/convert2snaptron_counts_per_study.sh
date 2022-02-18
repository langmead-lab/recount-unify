#!/usr/bin/env bash
set -exo pipefail
dir=$(dirname $0)
ddir="/datascope/recount03/rejoin_fix"
annot="G026"
ids="$ddir/recount3.sra.samples.tsv.cut"

study=$1

cwd=$(basename $PWD)
gene_sums=sra.gene_sums.${study}.${annot}.gz
if [[ $cwd == "tcga." || $cwd == "tcga.aa" ]]; then
    ids="$ddir/tcgav2.samples.tsv.cut"
    gene_sums=tcga.gene_sums.${study}.${annot}.gz
    if [[ $study == "NA" ]]; then
        ids="$ddir/tcgav2.samples.NA.tsv.cut"
    fi
fi
if [[ $cwd == "gtex." || $cwd == "gtex.aa" ]]; then
    ids="$ddir/gtexv2.samples.tsv.cut"
    gene_sums=gtex.gene_sums.${study}.${annot}.gz
fi
    
pushd $study

set +o pipefail
pcat $gene_sums | tail -n+3 | head -1 | cut -f 2- | tr '\t' '\n' > ${gene_sums}.samples
set -o pipefail

if [[ $cwd == "gtex." || $cwd == "gtex.aa" ]]; then
    /usr/bin/time -v pcat $gene_sums | tail -n+4 | python3 $dir/convert_gene_counts2snaptron.py ${gene_sums}.samples <(fgrep $'\t'"$study"$'\t' $ids | cut -f 1,3-4) | LC_ALL=C sort -k1,1 | cut -f2- > ${annot}.snaptron_full
else
    /usr/bin/time -v pcat $gene_sums | tail -n+4 | python3 $dir/convert_gene_counts2snaptron.py ${gene_sums}.samples <(fgrep $'\t'"$study"$'\t' $ids | cut -f 1-3) | LC_ALL=C sort -k1,1 | cut -f2- > ${annot}.snaptron_full
fi

echo "DONE"
popd
