#!/usr/bin/env bash
set -exo pipefail

dir=$(dirname $0)
orgn="human"
ddir=$PWD/../
ROOT=/datascope/recount03

study=$1
src=$2
threads=$3

if [[ -z $threads ]]; then
    threads=30
fi

#GP=$ROOT/release/${orgn}/data_sources/${src}/gene_sums
BP=$ROOT/release/${orgn}/data_sources/${src}/base_sums

pushd $study

lo=${study: -2}

#3 samples + gene
#num_samples2check=3

g=${src}.gene_sums.${study}.G026.gz
set +o pipefail
#num_samples=$(pcat $g | head -3 | tail -n1 | cut -f 2- | tr $'\t' $'\n' | wc -l)
#samples2check=$(pcat $g | head -3 | tail -n1 | tr $'\t' $'\n' | fgrep -n "" | egrep -v ':gene_id$' | shuf | tail -n${num_samples2check})
samples2check=$(pcat $g | head -3 | tail -n1 | tr $'\t' $'\n' | fgrep -n "" | egrep -v ':gene_id$')
set -o pipefail

BPP=$BP/$lo/$study

if [[ -z $SKIP_MD ]]; then
echo -n "" > resum.md.jobs
for idx_name in $samples2check ; do
    idx=$(echo $idx_name | cut -d':' -f1)
    sample=$(echo $idx_name | cut -d':' -f2)
            
    slo=${sample: -2}
    bw=$BPP/$slo/${src}.base_sums.${study}_${sample}.ALL.bw

    echo "$dir/megadepth120rc $bw --annotation $ddir/disjoint2exons2genes.G026.sorted.cut.bed --no-annotation-stdout --prefix $sample" >> resum.md.jobs
done
/usr/bin/time -v parallel -j${threads} < resum.md.jobs > resum.md.jobs.run 2>&1
fi

 
for idx_name in $samples2check ; do
    idx=$(echo $idx_name | cut -d':' -f1)
    sample=$(echo $idx_name | cut -d':' -f2)
            
    slo=${sample: -2}

    check_count=$(cut -f 4 ${sample}.annotation.tsv | perl -ne 'chomp; $s+=$_; END { print "$s\n"; }')

    new_count=$(pcat ${src}.gene_sums.${study}.G026.gz | tail -n+4 | cut -f $idx | perl -ne 'chomp; $s+=$_; END { print "$s\n"; }')

    if [[ $check_count -ne $new_count ]]; then
        echo "BAD_COUNT $sample $idx $study $src"
        exit -1
    else
        echo "GOOD_COUNT $sample $idx $study $src"
    fi 
done
echo "ALL_"$num_samples2check"_COUNT_CHECKS_GOOD $study $src"
popd
