#!/usr/bin/env bash
set -exo pipefail
#sum the bamcount (or megadepth) quantification of an additional split (non-overlapping) exon annotation (e.g. gencode v43)
#per-sample, at the gene level (summing split exons into whole gene sums)
#assumes both input files are the same order (so they gene_ids line up with the split-exon sums)
OUTPUT_PATH="additional_gene_sums"
#should be list of gene_id's matching # and order of additionall split-exon annotation BED used in Pump pipeline (in bamcount/megadepth), no header, single column (line can have multiple gene_ids separated by a "," for a single split exon)
#e.g. /shared-data/research/genomics/datasets/recount3/prep4refresh/annotation/human/split_exonsv43.bed.gene_ids
gene_idsF=$1
#zst'd compressed file of 4 columns (chr,start,end,bp_sum), no header
#e.g. DRR001175!DRP000425!hg38!sra.all1.annotation.tsv.zst
split_exons_sumsF=$2
output_path=$3
if [[ -z $output_path ]]; then
    output_path=$OUTPUT_PATH
fi
sample=$(basename $split_exons_sumsF)
sample=$(echo "$sample" | cut -d'!' -f 1)

if [[ ! -d $output_path ]]; then
    mkdir -p $output_path
fi
/usr/bin/time -v paste $gene_idsF <(zstd -cd $split_exons_sumsF | cut -f 4) | perl -ne 'chomp; $f=$_; ($g,$b)=split(/\t/,$f,-1); @g=split(/,/,$g,-1); for $g0 (@g) { $h{$g0}+=$b; } END { for $g (sort { $a cmp $b } keys %h) { print "$g\t".$h{$g}."\n"; }}' > $output_path/${sample}.summed
