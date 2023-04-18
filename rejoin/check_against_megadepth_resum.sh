#!/usr/bin/env bash
set -exo pipefail
#example run:
#REF_DIR=/container-mounts/recount/ref/hg38_unify /shared-data/research/genomics/software/monorail/recount-unify/rejoin/check_against_megadepth_resum.sh sra SRP141364 /container-mounts/recount/output/SRP141364_20220209_203246 `pwd` 3

dir=$(dirname $0)

#e.g. "sra"
src=$1
study=$2
#this will be the study-specific directory under which are all the sample
#sudirs and their base sums files
base_sums_parent_dir=$3
#e.g. `pwd`
root=$4
#optional
num_samples2check=$5
threads=$6
#G026, M023, whatever rat is
annot=$7
if [[ -z $annot ]];
    then annot="G026"
fi

pushd $root

if [[ -z $threads ]]; then
    threads=20
fi

lo=${study: -2}

#e.g. bw='/path/to/SRR7047911!SRP141364!hg38!sra.all.bw' 
g=$root/gene_sums_per_study/$lo/$study/${src}.gene_sums.${study}.${annot}.gz
set +o pipefail
pcat $g | head -3 | tail -n1 | tr $'\t' $'\n' | fgrep -n "" | egrep -v ':gene_id$' > ${study}.sample_idxs
set -o pipefail

if [[ -z $num_samples2check ]]; then
   num_samples2check=$(cat ${study}.sample_idxs | wc -l)
fi

set +o pipefail
bws=$(find $base_sums_parent_dir -name "*!${study}!*.all.bw" | head -${num_samples2check})
set -o pipefail

echo -n "" > samples_checked.txt
echo -n "" > resum.md.${study}.jobs
for bw in $bws; do
    bwb=$(basename $bw)
    sample=$(echo "$bwb" | cut -d'!' -f 1)
    slo=${sample: -2}
    echo "/usr/bin/time -v $dir/megadepth120 $bw --annotation $REF_DIR/disjoint2exons2genes.${annot}.sorted.cut.bed --no-annotation-stdout --prefix $sample > ${sample}.md_run 2>&1" >> resum.md.${study}.jobs
    echo "$sample" >> samples_checked.txt
done
/usr/bin/time -v parallel -j${threads} < resum.md.${study}.jobs > resum.md.${study}.jobs.run 2>&1

BAD=
for sample in `cat samples_checked.txt` ; do
    idx=$(egrep -e ':'"$sample"'$' ${study}.sample_idxs | cut -d':' -f1)
            
    slo=${sample: -2}

    check_count=$(cut -f 4 ${sample}.annotation.tsv | perl -ne 'chomp; $s+=$_; END { print "$s\n"; }')

    new_count=$(pcat $g | tail -n+4 | cut -f $idx | perl -ne 'chomp; $s+=$_; END { print "$s\n"; }')

    if [[ $check_count -ne $new_count ]]; then
        echo "BAD_COUNT $sample $idx $study $src"
        BAD=1
    else
        echo "GOOD_COUNT $sample $idx $study $src"
    fi 
done
if [[ -n $BAD ]]; then
        exit -1
fi
echo "ALL_"$num_samples2check"_COUNT_CHECKS_GOOD $study $src"
popd
