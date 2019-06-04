#!/usr/bin/env bash

#e.g. g26_g29_f6_r109.annotations.sorted.2019-02-18.gtf
joined_annotation_gtf=$1
#e.g. g26_g29_f6_r109.disjoint-exons.sorted.2019-02-18.bed
disjoint_exon_bed=$2

#create a bed file of unique exons from unioned annotated GTF
fgrep '	exon	' $joined_annotation_gtf | cut -f 1,4-7,9 | perl -ne 'BEGIN { %h=("G026"=>0,"G029"=>2, "R109"=>3, "F006"=>4); } chomp; $f=$_; @f=split(/\t/,$f); $f[5]=~s/^.*gene_id "([^"]+)".*$/$1/; @f1=split(/\./,$f[5]); $t=pop(@f1); $tid=$h{$t}; $g=pop(@f); $f[3]=$g."\t$tid"; $f[1]--; print "".join("\t",@f)."\n";' | sort -k1,1 -k2,2n -k3,3n -k5,5n | uniq  > ${joined_annotation_gtf}.exons.bed

#intersect the disjoint exons (a) with the union of all annotations' exons (b), must match strand (-s)
bedtools intersect -sorted -s -wao -a $disjoint_exon_bed -b ${joined_annotation_gtf}.exons.bed > disjoin.2019-02-18.gtf_exons.intersected.tsv

#create unique set of disjoint exons mapping to a list of >=1 genes
cat disjoin.2019-02-18.gtf_exons.intersected.tsv | cut -f 1-6,10  | perl -ne 'chomp; @f=split(/\t/,$_); $g=pop(@f); $f=join("\t",@f); if($f ne $pf) { if($pf) { print "$pf\t$pg\n"; } $pf=$f; $pg=$g; next; } $pg.=";$g"; END { if($pf) { print "$pf\t$pg\n"; } }' | gzip > disjoin.2019-02-18.gtf_exons.intersected.genes.tsv.gz

#create all exons (not disjoint) across all 4 annotations
#fgrep '	exon	' $joined_annotation_gtf | cut -f1,4-7,9 | perl -ne 'BEGIN { %h=("G026"=>0,"G029"=>2, "R109"=>3, "F006"=>4); } chomp; @f=split(/\t/,$_); $f[5]=~s/^.*gene_id "([^"]+)".*$/$1/; @f1=split(/\./,$f[5]); $t=pop(@f1); $tid=$h{$t}; $g=pop(@f); $f[3]=$g."\t$tid"; $f[1]--; print "".join("\t",@f)."\n";'  | sort -k1,1 -k2,2n -k3,3n -k5,5n > ${joined_annotation_gtf}.exons.sorted.bed

#get maximum extent of genes by their exons
#work around multiple copies of MIR219A1
sort -s -k4,4 ${joined_annotation_gtf}.exons.bed | perl -ne 'chomp; ($c,$s,$e,$g,$n,$o)=split(/\t/,$_); if($h{$g} && $h{$g}->[2] ne $c) { $g="$c"."_$g"; } $h{$g}->[0]=$s if(!defined($h{$g}->[0]) || $s < $h{$g}->[0]); $h{$g}->[1]=$e if(!defined($h{$g}->[1]) || $e > $h{$g}->[1]); $h{$g}->[2]=$c; $h{$g}->[3]=$o; END { for $g (keys %h) { ($s,$e,$c,$o)=@{$h{$g}}; print "$c\t$s\t$e\t$g\t0\t$o\n"; }}' | sort -k1,1 -k2,2n -k3,3n -k4,4 > g26_g29_f6_r109.annotations.sorted.2019-02-18.gtf.exons2genes.sorted.bed3

bedtools intersect -sorted -s -wao -a $disjoint_exon_bed -b g26_g29_f6_r109.annotations.sorted.2019-02-18.gtf.exons2genes.sorted.bed3 > ${joined_annotation_gtf}.disjoint2exons2genes.bed
