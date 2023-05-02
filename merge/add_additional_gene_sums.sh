#!/usr/bin/env
set -exo pipefail
sdir=$(dirname $0)

#/path/to/split_exonsv43.bed.gene_ids
gene_idsF=$1
#/path/to/exonsv43.bed.coords
exon_idsF=$2
#e.g. hg38
ref=$3
#e.g. gene_sums_per_study/<study_lo>/<study>/sra.gene_sums.DRP000425.G043.gz
output_file=$4

#optional
threads=$5

fn=$(basename $output_file)
#e.g. sra
src=$(echo "$fn" | cut -d'.' -f1)
#e.g. DRP000425
study=$(echo "$fn" | cut -d'.' -f3)
#e.g. G043
annotation=$(echo "$fn" | cut -d'.' -f4)

if [[ -z $threads ]]; then
    threads=8
fi

#assumes we're running as part of the Unifier
echo -n "" > additional_gene_sums.files
echo -n "" > additional_exon_sums.files
for d in `fgrep "/" links/attempts.moved.complete | fgrep "/${study}/"`; do
    #e.g. /container-mounts/input10/SRP277410/79/SRR12449879/SRR12449879_in3_att0
    sample=$(basename $(dirname $d))
    #DRR001175!DRP000425!hg38!sra.all1.annotation.tsv.zst
    #megadepth output version
    #echo "$d/${sample}!${study}!${ref}!${src}.all2.annotation.tsv.zst" >> additional_gene_sums.files
    #echo "$d/${sample}!${study}!${ref}!${src}.all1.annotation.tsv.zst" >> additional_exon_sums.files
    #bamcount output version
    echo "$d/${sample}!${study}!${ref}!${src}.all.exon_bw_count2.zst" >> additional_gene_sums.files
    echo "$d/${sample}!${study}!${ref}!${src}.all.exon_bw_count1.zst" >> additional_exon_sums.files
done
    
###Split Gene Sums:
#setup output file headers
output_file0=$(echo "$output_file" | sed 's#.gz$##')
echo "##annotation=$annotation" > $output_file0
date=$(date "+%Y%m%d %H:%m:%S")
#echo "##date.generated=2023-04-21 17:12:12.436322" >> $output_file0
echo "##date.generated=$date" >> $output_file0
sample_header="gene_id"
outputs=""
echo -n "" > additional_gene_sums.jobs
for f in `cat additional_gene_sums.files`; do
    sample=$(basename $f)
    sample=$(echo "$sample" | cut -d'!' -f 1)
    echo "/usr/bin/time -v /bin/bash -x $sdir/add_additional_gene_sums_per_sample.sh $gene_idsF $f additional_gene_sums/${sample}.summed > additional_sums_runs/${sample}.run 2>&1" >> additional_gene_sums.jobs
    sample_header="${sample_header}	${sample}"
    outputs="$outputs additional_gene_sums/${sample}.summed"
done
echo "$sample_header" >> $output_file0
mkdir -p additional_gene_sums
mkdir -p additional_sums_runs
/usr/bin/time -v parallel -j${threads} < additional_gene_sums.jobs > additional_gene_sums.jobs.run${threads} 2>&1
#now paste the results together
paste $gene_idsF $outputs >> $output_file0
pigz -f --fast -p8 $output_file0


###Full exon sums (no need to do summing, just pasting):
output_file0=$(echo "$output_file" | sed 's#.gz$##' | sed 's#gene#exon#g')
echo "##annotation=$annotation" > $output_file0
echo "##date.generated=$date" >> $output_file0
sample_header="chromosome|start_1base|end_1base|strand"
outputs=""
echo -n "" > exon_uzstd.jobs
mkdir -p additional_exon_sums
for f in `cat additional_exon_sums.files`; do
    sample=$(basename $f)
    sample=$(echo "$sample" | cut -d'!' -f 1)
    echo "/usr/bin/time -v zstd -cd $f 2> additional_sums_runs/${sample}.exon.run | cut -f 4 > additional_exon_sums/${sample}.summed" >> exon_uzstd.jobs 
    sample_header="${sample_header}	${sample}"
    f0=$(echo "$f" | sed 's#.zst$##')
    outputs="$outputs additional_exon_sums/${sample}.summed"
done
echo "$sample_header" >> $output_file0
/usr/bin/time -v parallel -j${threads} < exon_uzstd.jobs > exon_uzstd.jobs.run${threads} 2>&1
paste $exon_idsF $outputs >> $output_file0
pigz -f --fast -p8 $output_file0
