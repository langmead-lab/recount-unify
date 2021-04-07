#!/usr/bin/env bash
set -exo pipefail

dir=$(dirname $0)

#e.g. srav3_human5
tranche=$1
#disjoint2exons.bed
#exon_mapping_file=$2
#no more than (for sorting/compression)
threads=$2

#input="all.exon_bw_count.pasted.gz"
output_old="all.exon_counts.rejoined.tsv.old_ordering.gz"
output="all.exon_counts.rejoined.tsv.gz"

pushd $tranche
#pcat is pigz 2-thread version of zcat
#head causes pipefail to crash the script
#set +o pipefail
#num_samples=$(pcat $input | head -1 | cut -f 7- | tr \\t \\n | wc -l)
#set -o pipefail
#echo "now rejoining from $dir $exon_mapping_file $input $num_samples"
#$dir/rejoin -a ${exon_mapping_file} -d <(pcat ${input}) -s ${num_samples} -p exon -h
export LC_ALL=C
#sort -t'	' -k2,2 -k3,3n -k4,4n -k6,6 --stable --parallel=${threads} exon.counts > exon.counts.sorted
pcat $output_old > all.exon_counts.rejoined.old_ordering.tsv
sort -T/home/cwilks3/langmead/temp -t'	' -k2,2 -k3,3n -k4,4n -k6,6 --stable --parallel=${threads} all.exon_counts.rejoined.old_ordering.tsv > exon.counts.sorted
#rm -f exon.counts exon.intron_counts
cut -f 1-6 exon.counts.sorted > ${output}.coords
#cat exon.counts.sorted | pigz --fast -p ${threads} > ${output}
pigz --fast -p ${threads} exon.counts.sorted
mv exon.counts.sorted.gz $output
#rm -f exon.counts.sorted
rm -f all.exon_counts.rejoined.old_ordering.tsv
popd
