#!/usr/bin/env
set -exo pipefail
#ulimit -n 25000
dir=$(dirname $0)
orgn="human"
ddir=$PWD/../
ROOT=/datascope/recount03
rejoin=$ROOT/recount-unify/rejoin/rejoin.postfix
bed=$ddir/disjoint2exons2genes.fix.bed.disjoint_exons.sorted
sanity_check=$ddir/good_genes_sanity_check.txt
#bed=$dir/recount3.exons.bed

study=$1
src=$2

GP=$ROOT/release/${orgn}/data_sources/${src}/gene_sums
BP=$ROOT/release/${orgn}/data_sources/${src}/base_sums

mkdir -p $study
pushd $study

bed_short=$(basename $bed)

lo=${study: -2}

for annot in G026 G029 R109 F006; do
set +o pipefail
    pcat $GP/$lo/$study/${src}.gene_sums.${study}.${annot}.gz | head -3 | tail -n1 | tr $'\t' $'\n' | tail -n+2 > ${annot}.samples
set -o pipefail
done

num_samples=$(cat G026.samples | wc -l)
BPP=$BP/$lo/$study
if [[ -z $SKIP_MD ]]; then
    echo -n "" > list_of_output.files
    echo -n "" > header.1
    echo -n "" > md.jobs
    n=0
    for sample in `cat G026.samples`; do
        n=$((n + 1))
        slo=${sample: -2}
        #bw=$BPP/$slo/$sample/${src}.base_sums.${study}_${sample}.ALL.bw
        bw=$BPP/$slo/${src}.base_sums.${study}_${sample}.ALL.bw
        if [[ $src == "gtex" ]]; then
            slo=${sample: -4}
            slo=$(echo "$slo" | sed 's#..$##')
            bw=$BPP/$slo/${src}.base_sums.${study}_${sample}.ALL.bw
        fi
        if [[ $src == "tcga" ]]; then
            slo=${slo^^}
            bw=$BPP/$slo/${src}.base_sums.${study}_${sample}.ALL.bw
        fi
        echo "/usr/bin/time -v /bin/bash -x $dir/bigwig_coverage.sh $bw $bed ${sample}.md > ${sample}.md.run 2>&1" >> md.jobs
        echo "${sample}.md.notsummed " >> list_of_output.files
        echo -n "	$n" >> header.1
    done
    echo "" >> header.1
    /usr/bin/time -v parallel -j20 < md.jobs > md.jobs.run 2>&1
fi

set +o pipefail
pcat $GP/${src}.gene_sums.${study}.${annot}.gz | head -3 > ${src}.gene_sums.${study}.${annot}
set -o pipefail
list=$(cat list_of_output.files)
cat <(echo -n $'gene\tstart\tend\tname\tscore\tstrand') <(cat header.1) > ${study}.genes.pasted.tsv
#alternate approach to pasting to get around number of files open/too long command line
#from https://unix.stackexchange.com/questions/205642/combining-large-amount-of-files
#echo -n "" > ${study}.genes.pasted.tsv_
#set +e
#for f in $list; do cat ${study}.genes.pasted.tsv_ | paste - $f >temp; cp temp ${study}.genes.pasted.tsv_; done; rm temp
#cat ${study}.genes.pasted.tsv_ >> ${study}.genes.pasted.tsv ; rm ${study}.genes.pasted.tsv_
#set -e
paste $bed $list >> ${study}.genes.pasted.tsv

$rejoin -a $ddir/disjoint2exons2genes.fix.sorted.bed -d ${study}.genes.pasted.tsv -s $num_samples -p gene -h

for annot in G026 G029 R109 F006; do
    fgrep ".${annot}	" gene.counts | cut -f 1,7- > gene.counts.${annot}
    #then do update python3 script of original gene counts
    pcat $GP/$lo/$study/${src}.gene_sums.${study}.${annot}.gz | python3 $dir/update_counts.py <(cat gene.counts.${annot} | sed 's#\.'$annot'##') <(fgrep ".${annot}" $sanity_check | sed 's#\.'$annot'##' | sed 's#\t#|#g') 2> ${src}.gene_sums.${study}.${annot}.gz.update | pigz --fast -p2 > ${src}.gene_sums.${study}.${annot}.gz
    set +e
    diff <(pcat ${src}.gene_sums.${study}.${annot}.gz) <(pcat $GP/$lo/$study/${src}.gene_sums.${study}.${annot}.gz) > ${annot}.diff
    set -e
    fgrep "<" ${annot}.diff | wc -l > ${annot}.diff.check
    fgrep ">" ${annot}.diff | wc -l >> ${annot}.diff.check
done

popd
