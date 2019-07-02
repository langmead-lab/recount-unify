#!/usr/bin/env bash

#e.g. g26_g29_f6_r109.annotations.sorted.2019-02-18.gtf
joined_annotation_gtf=$1
#e.g. g26_g29_f6_r109.disjoint-exons.sorted.2019-02-18.bed
disjoint_exon_bed=$2

#create a bed file of unique exons from unioned annotated GTF
fgrep '	exon	' $joined_annotation_gtf | cut -f 1,4-7,9 | perl -ne 'BEGIN { %h=("G026"=>0,"G029"=>2, "R109"=>3, "F006"=>4); } chomp; $f=$_; @f=split(/\t/,$f); $f[5]=~s/^.*gene_id "([^"]+)".*$/$1/; @f1=split(/\./,$f[5]); $t=pop(@f1); $tid=$h{$t}; $g=pop(@f); $f[3]=$g."\t$tid"; $f[1]--; print "".join("\t",@f)."\n";' | sort -k1,1 -k2,2n -k3,3n -k5,5n | uniq  > ${joined_annotation_gtf}.exons.bed

#intersect the disjoint exons (a) with the union of all annotations' exons (b), must match strand (-s)
bedtools intersect -sorted -s -wao -a $disjoint_exon_bed -b ${joined_annotation_gtf}.exons.bed > disjoin.2019-02-18.gtf_exons.intersected.tsv

#UPDATE attempt at fixing extra disjoint exons being included which aren't fully contained within an originally annotated exon
#ultimately this is still in question, need a more robust mapping from Leo
#create unique set of disjoint exons mapping to a list of >=1 genes
cat gencode.v25.annotation.gtf.exons.bed2disjoin_exons_recount2_intersected  | perl -ne 'chomp; ($c1,$s1,$e1,$g1,$bp1,$o1,$c2,$s2,$e2,$g2,$j1,$o2)=split(/\t/,$_); next if(!($s2<=$s1 && $e2>=$e1)); $g=$g2; $f=join("\t",($c1,$s1,$e1,$g1,$bp1,$o1)); if($f ne $pf) { if($pf) { print "$pf\t$pg\n"; } $pf=$f; $pg=$g; next; } $pg.=";$g"; END { if($pf) { print "$pf\t$pg\n"; } }' > ${joined_annotation_gtf}.exons.intersected.genes

#create all exons (not disjoint) across all 4 annotations
#fgrep '	exon	' $joined_annotation_gtf | cut -f1,4-7,9 | perl -ne 'BEGIN { %h=("G026"=>0,"G029"=>2, "R109"=>3, "F006"=>4); } chomp; @f=split(/\t/,$_); $f[5]=~s/^.*gene_id "([^"]+)".*$/$1/; @f1=split(/\./,$f[5]); $t=pop(@f1); $tid=$h{$t}; $g=pop(@f); $f[3]=$g."\t$tid"; $f[1]--; print "".join("\t",@f)."\n";'  | sort -k1,1 -k2,2n -k3,3n -k5,5n > ${joined_annotation_gtf}.exons.sorted.bed

#get maximum extent of genes by their exons
#work around multiple copies of MIR219A1
sort -s -k4,4 ${joined_annotation_gtf}.exons.bed | perl -ne 'chomp; ($c,$s,$e,$g,$n,$o)=split(/\t/,$_); if($h{$g} && $h{$g}->[2] ne $c) { $g="$c"."_$g"; } $h{$g}->[0]=$s if(!defined($h{$g}->[0]) || $s < $h{$g}->[0]); $h{$g}->[1]=$e if(!defined($h{$g}->[1]) || $e > $h{$g}->[1]); $h{$g}->[2]=$c; $h{$g}->[3]=$o; END { for $g (keys %h) { ($s,$e,$c,$o)=@{$h{$g}}; print "$c\t$s\t$e\t$g\t0\t$o\n"; }}' | sort -k1,1 -k2,2n -k3,3n -k4,4 > g26_g29_f6_r109.annotations.sorted.2019-02-18.gtf.exons2genes.sorted.bed3

#TODO test the updated approach with the 4 annotations
#now do gene mapping from disjoint exons
#create a bed file of unique exons from unioned annotated GTF
#egrep -e '	gene	' $joined_annotation_gtf | cut -f 1,4,5,7,9 | cut -d'"' -f 1,2 | sed 's/gene_id "//' | perl -ne 'chomp; @f=split(/\t/,$_); $gid=pop(@f); $f[1]--; $f[2].="\t$gid\t0"; print "".join("\t",@f)."\t2\n";' | sort -k1,1 -k2,2n -k3,3n > ${joined_annotation_gtf}.genes
fgrep '	gene	' $joined_annotation_gtf | cut -f 1,4-7,9 | perl -ne 'BEGIN { %h=("G026"=>0,"G029"=>2, "R109"=>3, "F006"=>4); } chomp; $f=$_; @f=split(/\t/,$f); $f[5]=~s/^.*gene_id "([^"]+)".*$/$1/; @f1=split(/\./,$f[5]); $t=pop(@f1); $tid=$h{$t}; $g=pop(@f); $f[3]=$g."\t$tid"; $f[1]--; print "".join("\t",@f)."\n";' | sort -k1,1 -k2,2n -k3,3n -k5,5n | uniq  > ${joined_annotation_gtf}.genes.bed

cat ${joined_annotation_gtf}.genes.bed ${joined_annotation_gtf}.exons.intersected.genes | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); $gid=$f[3]; if($f[6] eq "2") { $h{$gid}=$f; next; } @gids=split(/;/,$f[6]); pop(@f); $f=join("\t",@f); %seen=(); for $gid (@gids) { next if($seen{$gid}); $ginfo=$h{$gid}; @f2=split(/\t/,$ginfo); pop(@f2); $ginfo=join("\t",@f2); print "$f\t$ginfo\t".$f[4]."\n"; $seen{$gid}=1; }' > ${joined_annotation_gtf}.exons.disjoint2exons2genes.bed
