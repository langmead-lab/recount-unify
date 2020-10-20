#!/usr/bin/env bash
set -euxo pipefail

dir=$(dirname $0)

#MM file (matrix market formatted jx file w/ 3 header lines), gzipped
#e.g. KIDNEY.all.MM.gz
mmf=$1
#e.g. KIDNEY.all.RR.gz
#RR file (jx row metadata w/ 1 header line), gzipped
rrf=$2
#list of rail_ids for all samples in the MM file, not gzipped
#e.g. KIDNEY.unique.sj.merged.motifs.annotated.unique.sids
#KIDNEY.all.sj.merged.motifs.annotated.sids
#sids=$3
#tcga, gtex, sra
comp=$3

#works for gtex, not for tcga
sids=$(echo $mmf | perl -ne 'chomp; $f=$_; @f=split(/\./,$f); $d=shift(@f); $t=shift(@f); $r="$d.$t.sj.merged.motifs.annotated.$t.sids"; if($t eq "all") { $r=~s/.annotated.all.sids/.annotated.sids/; } print "$r\n";')
if [[ "$comp" != "gtex" ]]; then
    sids=$(echo $mmf | perl -ne 'chomp; $f=$_; @f=split(/\./,$f); $d=shift(@f); $t=shift(@f); $r="$d.unique.sj.merged.motifs.annotated.unique.sids"; print "$r\n";')
fi

{ pcat $mmf ||:; } | head -3 | tail -n1 | tr \\t \\n > ${mmf}.check

num_rows=$(head -1 ${mmf}.check)
num_cols=$(head -2 ${mmf}.check | tail -n1)
expected_num_vals=$(head -3 ${mmf}.check | tail -n1)

#skip header lines
pcat $mmf | tail -n+4 | cut -f 2 | ${dir}/count_n_check -m $num_cols >> ${mmf}.check

num_max_cols=$(head -4 ${mmf}.check | tail -n1)
actual_num_vals=$(head -5 ${mmf}.check | tail -n1)

pcat $rrf | tail -n+2 | wc -l >> ${mmf}.check

num_jxs=$(head -6 ${mmf}.check | tail -n1)

cat $sids | wc -l >> ${mmf}.check

num_sids=$(tail -n1 ${mmf}.check)

#double print output, lazy way of capturing it in both ways depending on how the script is run
if [[ $num_rows -ne $num_jxs ]]; then
    echo "ERROR	$mmf	num_rows != num_jxs: $num_rows != $num_jxs" >> ${mmf}.check
    echo "ERROR	$mmf	num_rows != num_jxs: $num_rows != $num_jxs"
    exit -1
fi

if [[ $num_cols -ne $num_sids ]]; then
    echo "ERROR	$mmf	num_cols != num_sids: $num_cols != $num_sids" >> ${mmf}.check
    echo "ERROR	$mmf	num_cols != num_sids: $num_cols != $num_sids"
    exit -1
fi

if [[ $expected_num_vals -ne $actual_num_vals ]]; then
    echo "ERROR	$mmf	expected_num_vals != actual_num_vals: $num_rows != $actual_num_vals" >> ${mmf}.check
    echo "ERROR	$mmf	expected_num_vals != actual_num_vals: $num_rows != $actual_num_vals"
    exit -1
fi

#max_col_id_in_file=$(pcat $mmf | fgrep -m2 "	$num_cols	" | wc -l)

#if [[ $max_col_id_in_file -ne 2 ]]; then
if [[ $num_max_cols -lt 2 ]]; then
    echo "ERROR	$mmf	max column ID not found in MM file: $num_max_cols" >> ${mmf}.check
    echo "ERROR	$mmf	max column ID not found in MM file: $num_max_cols"
    exit -1
fi

echo "$mmf	OK" >> ${mmf}.check
echo "$mmf	OK"
