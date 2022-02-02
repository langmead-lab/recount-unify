#!/usr/bin/env
set -exo pipefail
dir=$(dirname $0)
ROOT=
GP=$ROOT/release/human/data_sources/sra/gene_sums
BP=$ROOT/release/human/data_sources/sra/base_sums
rejoin=$ROOT/recount-unify/rejoin/rejoin.postfix
#bed=$dir/disjoint2exons2genes.fix.bed.disjoint_exons.sorted
bed=$dir/recount3.exons.bed
orgn=human

study=$1
src=$2
annot=$3

if [[ -z $annot ]]; then
    annot="G026"
fi

mkdir -p $study
pushd $study

bed_short=$(basename $bed)

lo=${study: -2}

set +o pipefail
pcat $GP/$lo/$study/${src}.gene_sums.${study}.${annot}.gz | head -3 | tail -n1 | tr $'\t' $'\n' | tail -n+2 > samples
set -o pipefail
num_samples=$(cat samples | wc -l)
BPP=$BP/$lo/$study
if [[ -z $SKIP_MD ]]; then
    echo -n "" > list_of_output.files
    echo -n "" > header.1
    echo -n "" > md.jobs
    n=0
    for sample in `cat samples`; do
        n=$((n + 1))
        slo=${sample: -2}
        bw=$BPP/$slo/${src}.base_sums.${study}_${sample}.ALL.bw
        echo "/usr/bin/time -v /bin/bash -x $dir/bigwig_coverage.sh $bw $bed ${sample}.md > ${sample}.md.run 2>&1" >> md.jobs
        echo "${sample}.md.notsummed " >> list_of_output.files
        echo -n "	$n" >> header.1
    done
fi

/usr/bin/time -v parallel -j20 < md.jobs > md.jobs.run 2>&1

set +o pipefail
pcat $GP/${src}.gene_sums.${study}.${annot}.gz | head -3 > ${src}.gene_sums.${study}.${annot}
set -o pipefail
list=$(cat list_of_output.files)
cat <(echo -n $'gene\tstart\tend\tname\tscore\tstrand') <(cat header.1) > ${study}.genes.pasted.tsv
paste $bed $list >> ${study}.genes.pasted.tsv

$rejoin -a $dir/disjoint2exons2genes.fix.sorted.bed -d ${study}.genes.pasted.tsv -s $num_samples -p gene -h

popd
