#!/usr/bin/env bash
set -exo pipefail
dir=$(dirname $0)
#given a BED file (6 columns) and a BigWig file
#calculate the 1) gene and 2) exon coverage for the BED file's regions (assumed to be exons) over the base-level coverage in the BigWig
#this assumes:
#1) BED file is sorted by coordinate, e.g. sort -k1,1 -k2,2n
#2) spliced regions for a specific gene are contiguous in BED file (e.g. for a gene, it's exonic regions are not interrupted by another gene's exon)
#3) spliced genes are marked about by their exon regions containing in the BED file's name column (4th column) >=1 characters and the suffix: ":1"

#this script will first calculate the raw exon coverage using Megadepth
#then this script will sum the coverage columns between exons to get gene level coverages for the sample

#MD=$dir/megadepth113
#faster than 113 and correct
MD=$dir/megadepth120rc
#genes=$dir/disjoint2exons2genes.fix.sorted.bed.exons.bed.genes
#bed=$dir/disjoint2exons2genes.fix.sorted.bed.exons.bed
#bed=$dir/disjoint2exons2genes.fix.bed.disjoint_exons.sorted

bw=$1
bed=$2
prefix=$3
sum=1

$MD $bw --annotation $bed --no-annotation-stdout --sums-only --prefix $prefix
if [[ -n $sum ]]; then
    paste $genes ${prefix}.annotation.tsv > ${prefix}.pasted
    cat ${prefix}.pasted | perl -ne 'chomp; $f=$_; ($n,$v)=split(/\t/,$f,-1); $v=~s/\.\d+$//; if(defined($h{$n})) { $h{$n}+=$v; next; } $h{$n}=$v; END { for $g (keys %h) { print "$g\t".$h{$g}."\n"; }}' > ${prefix}.summed
    LC_ALL=C sort -k1,1 ${prefix}.summed > ${prefix}.summed.sorted
    cut -f 2 ${prefix}.summed.sorted > ${prefix}.summed
fi
cat ${prefix}.annotation.tsv | perl -ne 'chomp; $f=$_; $f=~s/\.\d+$//; print "$f\n";' > ${prefix}.notsummed
rm ${prefix}.annotation.tsv
