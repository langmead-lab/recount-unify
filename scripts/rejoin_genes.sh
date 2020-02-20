#!/usr/bin/bash

#paste disjoint exon sums file
input=$1
#file to write output to (gzipped)
output=$2
#disjoint exon2annotated[gene/exon] file
mapping=$3
#path to rejoin binary
rejoin_path=$4
#type: gene or exon
t=$5
#[optional] number of threads to use when pigz'ing output, default 4
threads=$6
if [[ -z $threads ]]; then
	threads=4
fi

#rejoin binary needs exact number of runs/samples
NUM_SAMPLES=`zcat ${input} | head -1 | cut -f 7- | tr \\\\t \\\\n | wc -l`
#actual call to rejoin binary, most of the time intensive stuff happens here
${rejoin_path} -a ${mapping} -d <(zcat ${input}) -s $NUM_SAMPLES -p $t -p ${output} -h
#stream header in along with sums and compress
cat <(echo -n "gene	chromosome	start	end	length	strand	") <(zcat ${input} | head -1 | cut -f 7-) ${output}.counts | pigz --fast -p $threads > ${output}
#dont need intermediate files
rm -f ${output}.counts ${output}.intron_counts
