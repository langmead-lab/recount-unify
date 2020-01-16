#!/usr/bin/bash

input=$1
output=$2
mapping=$3
rejoin_path=$4
t=$5

NUM_SAMPLES=`zcat ${input} | head -1 | cut -f 7- | tr \\\\t \\\\n | wc -l`
${rejoin_path} -a ${mapping} -d <(zcat ${input}) -s $NUM_SAMPLES -p $t -p ${output} -h
cat <(echo -n "gene	chromosome	start	end	length	strand	") <(zcat ${input} | head -1 | cut -f 7-) ${output}.counts | gzip > ${output}
rm -f ${output}.counts ${output}.intron_counts
