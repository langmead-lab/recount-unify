#!/usr/bin/env
set -exo pipefail
dir=$(dirname $0)
orgn="human"
ddir=$PWD/../
ROOT=/datascope/recount03
rejoin=$ROOT/recount-unify/rejoin/rejoin.postfix
bed=$ddir/disjoint2exons2genes.fix.bed.disjoint_exons.sorted
num_bed_lines=146725
sanity_check=$ddir/good_genes_sanity_check.txt
#bed=$dir/recount3.exons.bed
threads=20

study=$1
src=$2

GP=$ROOT/release/${orgn}/data_sources/${src}/gene_sums
BP=$ROOT/release/${orgn}/data_sources/${src}/base_sums

pushd $study

bed_short=$(basename $bed)
lo=${study: -2}

num_samples=$(cat G026.samples | wc -l)
BPP=$BP/$lo/$study
REDO_DIFFS=
if [[ ! -f pasted.outputs ]]; then
    wc -l *.md.notsummed  > line_checks
    num_good_files=$(fgrep "$num_bed_lines" line_checks | egrep -v -e " total$" | wc -l)
    if [[ $num_good_files -ne $num_samples ]]; then
        echo "not all samples have required number of lines, terminating early!"
        exit -1
    fi 
    MAX_FILES_PER_PASTE=500
    num_files=$(cat list_of_output.files | wc -l)
    cat list_of_output.files | perl -ne 'BEGIN { $max='$MAX_FILES_PER_PASTE'; } chomp; $f=$_; $i++; $files.="$f "; if($i % $max == 0) { print "paste $files > pasted.$i.files\n"; print STDERR "pasted.$i.files\n"; $files=""; } END { if(length($files) > 0) { $i++; print "paste $files > pasted.$i.files\n"; print STDERR "pasted.$i.files\n"; }}' > paste.jobs 2>pasted.outputs
    /usr/bin/time -v parallel -j$threads < paste.jobs > paste.jobs.run 2>&1
    #list=$(cat list_of_output.files)
    list=$(cat pasted.outputs)
    cat <(echo -n $'gene\tstart\tend\tname\tscore\tstrand') <(cat header.1) > ${study}.genes.pasted.tsv
    paste $bed $list >> ${study}.genes.pasted.tsv
    $rejoin -a $ddir/disjoint2exons2genes.fix.sorted.bed -d ${study}.genes.pasted.tsv -s $num_samples -p gene -h
    REDO_DIFFS=1
fi

for annot in G026 G029 R109 F006; do
    if [[ -n $REDO_DIFFS || ! -f ${annot}.diff ]]; then
        fgrep ".${annot}	" gene.counts | cut -f 1,7- > gene.counts.${annot}
        #then do update python3 script of original gene counts
        pcat $GP/$lo/$study/${src}.gene_sums.${study}.${annot}.gz | python3 $dir/update_counts.py <(cat gene.counts.${annot} | sed 's#\.'$annot'##') <(fgrep ".${annot}" $sanity_check | sed 's#\.'$annot'##' | sed 's#\t#|#g') 2> ${src}.gene_sums.${study}.${annot}.gz.update | pigz --fast -p2 > ${src}.gene_sums.${study}.${annot}.gz
        set +e
        diff <(pcat ${src}.gene_sums.${study}.${annot}.gz) <(pcat $GP/$lo/$study/${src}.gene_sums.${study}.${annot}.gz) > ${annot}.diff
        set -e
    fi
    set +o pipefail
    num_diff_lines=$(fgrep "<" ${annot}.diff | wc -l)
    set -o pipefail
    if [[ $num_diff_lines -gt 0 ]]; then
        expected_total_num_fields=$((num_samples * num_diff_lines)) 
        total_num_fields=$(fgrep "<" ${annot}.diff | cut -f2- | tr $'\t' $'\n' | wc -l)
        if [[ $total_num_fields -ne $expected_total_num_fields ]]; then
            echo "bad number of fields in new part of diff, terminating early!"
            exit -1
        fi
        echo "$num_diff_lines" > ${annot}.diff.check
        fgrep ">" ${annot}.diff | wc -l >> ${annot}.diff.check
    fi
done
#check that we get expected number of diffs and that they agree between the two files
if [[ -f G026.diff.check ]]; then
    cat G026.diff.check | tr $'\n' ':' | perl -ne 'chomp; ($n1,$n2)=split(/:/,$_,-1); die "bad diff count $n1 for G026\n" if($n1 != $n2 || $n1 > 1797);'
fi
if [[ -f G029.diff.check ]]; then
    cat G029.diff.check | tr $'\n' ':' | perl -ne 'chomp; ($n1,$n2)=split(/:/,$_,-1); die "bad diff count $n1 for G029\n" if($n1 != $n2 || $n1 > 1364);'
fi
#was 74
if [[ -f F006.diff.check ]]; then
    cat F006.diff.check | tr $'\n' ':' | perl -ne 'chomp; ($n1,$n2)=split(/:/,$_,-1); die "bad diff count $n1 for F006\n" if($n1 != $n2 || $n1 > 80);'
fi
if [[ -f R109.diff.check ]]; then
    cat R109.diff.check | tr $'\n' ':' | perl -ne 'chomp; ($n1,$n2)=split(/:/,$_,-1); die "bad diff count $n1 for R109\n" if($n1 != $n2 || $n1 > 2589);'
fi

popd
