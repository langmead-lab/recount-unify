#e.g. G026
src=$1
#e.g. G029
#target=$2
export LC_ALL=C

src_short=$src
for target in G026 G029 R109 F006; do

src=$src_short
#support self comparison
targeto=$target
if [[ $src == $target ]]; then
    target=${target}"b"
fi

fgrep "$src" ../disjoint2exons2genes.bed > disjoint2exons2genes.bed.${src}
fgrep "$targeto" ../disjoint2exons2genes.bed > disjoint2exons2genes.bed.${target}

src_short=$src
target_short=$target

src=disjoint2exons2genes.bed.${src}
sort -k10,10 -k1,1 -k2,2n -k3,3n -k6,6 ${src}  > ${src}.sorted2
cat ${src}.sorted2 |perl -ne 'chomp; $f=$_; @f=split(/\t/,$f,-1); ($c,$s,$e)=@f; $o=$f[5]; $gc=$f[6]; $gs=$f[7]; $ge=$f[8]; $go=$f[11]; $g=$f[9]; if($pg) { if($pg eq $g) { $pes.=";$s-$e"; $ne++; next; } print "$pg\t$ne\t$pc\t$ps\t$pe\t$po\t$pes\n"; } $ne=1; $pg=$g; $pc=$gc; $ps=$gs; $pe=$ge; $po=$go; $pes="$s-$e"; END { if($pg) { print "$pg\t$ne\t$pc\t$ps\t$pe\t$po\t$pes\n"; }}' > ${src}.exonstrings

target=disjoint2exons2genes.bed.${target}
sort -k10,10 -k1,1 -k2,2n -k3,3n -k6,6 ${target} > ${target}.sorted2
#cat ${target}.sorted2 |perl -ne 'chomp; $f=$_; @f=split(/\t/,$f,-1); ($c,$s,$e)=@f; $o=$f[5]; $gc=$f[6]; $gs=$f[7]; $ge=$f[8]; $go=$f[11]; $g=$f[9]; if($pg) { if($pg eq $g) { $pes.=";$s-$e"; next; } print "$pg\t$pc\t$ps\t$pe\t$po\t$pes\n"; } $pg=$g; $pc=$gc; $ps=$gs; $pe=$ge; $po=$go; $pes="$s-$e"; END { if($pg) { print "$pg\t$pc\t$ps\t$pe\t$po\t$pes\n"; }}' > ${target}.exonstrings
cat ${target}.sorted2 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f,-1); ($c,$s,$e)=@f; $o=$f[5]; $gc=$f[6]; $gs=$f[7]; $ge=$f[8]; $go=$f[11]; $g=$f[9]; if($pg) { if($pg eq $g) { $pes.=";$s-$e"; $ne++; next; } print "$pg\t$ne\t$pc\t$ps\t$pe\t$po\t$pes\n"; } $ne=1; $pg=$g; $pc=$gc; $ps=$gs; $pe=$ge; $po=$go; $pes="$s-$e"; END { if($pg) { print "$pg\t$ne\t$pc\t$ps\t$pe\t$po\t$pes\n"; }}' > ${target}.exonstrings

cut -f 3-6 ${src}.exonstrings > ${src}.exonstrings.coords

fgrep -f ${src}.exonstrings.coords ${target}.exonstrings | cut -f 3- | sort -k1,1 -k2,2n -k3,3n -k4,4 > ${target}.exonstrings.by${src_short}

cut -f 1-4 ${target}.exonstrings.by${src_short} > ${target}.exonstrings.by${src_short}.coords
fgrep -f ${target}.exonstrings.by${src_short}.coords ${src}.exonstrings | cut -f 3- | sort -k1,1 -k2,2n -k3,3n -k4,4 > ${src}.exonstrings.bycoords

#diff ${src}.exonstrings.bycoords ${target}.exonstrings.by${src_short} > ${src_short}.vs.${target_short}.diff
#fgrep "<" ${src_short}.vs.${target_short}.diff | cut -d' ' -f 2- | cut -f 1-4 > ${src_short}.vs.${target_short}.diff.${src_short}
#instead of simple diff, do hash check, still exact matching but more flexible
cat ${target}.exonstrings.by${src_short} | perl -ne 'BEGIN { open(IN,"<'${src}.exonstrings.bycoords'"); while($line=<IN>) { chomp($line); ($c,$s,$e,$o,$p)=split(/\t/,$line,-1); $k=join("|",($c,$s,$e,$o)); push(@{$h{$k}},$p); } close(IN); } chomp; $f=$_; ($c,$s,$e,$o,$p2)=split(/\t/,$f,-1); $k=join("|",($c,$s,$e,$o)); $p_=$h{$k}; if(defined($p_)) { for $p00 (@$p_) { if($p2 ne $p00) { $k=~s/\|/\t/g; print "$k\n"; }}}' > ${src_short}.vs.${target_short}.diff
ln -fs ${src_short}.vs.${target_short}.diff ${src_short}.vs.${target_short}.diff.${src_short}

fgrep -f ${src_short}.vs.${target_short}.diff.${src_short} ${src}.exonstrings | sort -u | sort -k3,3 -k4,4n -k5,5n > ${src}.exonstrings.${src_short}.vs.${target_short}
fgrep -f ${src_short}.vs.${target_short}.diff.${src_short} ${target}.exonstrings | sort -u | sort -k3,3 -k4,4n -k5,5n > ${target}.exonstrings.${src_short}.vs.${target_short}

#cat ${target}.exonstrings.${src_short}.vs.${target_short} | perl -ne 'BEGIN { open(IN,"<'${src}.exonstrings.${src_short}.vs.${target_short}'"); while($line=<IN>) { chomp($line); ($g,$n,$c,$s,$e,$o,$exons)=split(/\t/,$line); @f2=split(/;/,$exons); $k=join("|",($c,$s,$e,$o)); for $ex (@f2) { $h{$k}->{$g}->{$ex}=1; }} close(IN); } chomp; $f=$_; ($g,$n,$c,$s,$e,$o,$exons)=split(/\t/,$f); @f2=split(/;/,$exons); $k=join("|",($c,$s,$e,$o)); $n_genes_contained=0; @keys=keys %{$h{$k}}; $n_othergenes=scalar(@keys); for $othergene (@keys) { if($othergene eq $g) { $n_othergenes--; next; } $n3=scalar(keys %{$h{$k}->{$othergene}}); next if($n3 > $n); $n2=0; for $exon (@f2) { if(defined($h{$k}->{$othergene}->{$exon})) { $n2++; } }  if($n2 == $n3) { $n_genes_contained++; } } if($n_othergenes == $n_genes_contained) { print "contained\t$f\n"; }' > ${target}.exonstrings.${src_short}.vs.${target_short}.contained
cat ${src}.exonstrings.${src_short}.vs.${target_short} | perl -ne 'BEGIN { open(IN,"<'${target}.exonstrings.${src_short}.vs.${target_short}'"); while($line=<IN>) { chomp($line); ($g,$n,$c,$s,$e,$o,$exons)=split(/\t/,$line); @f2=split(/;/,$exons); $k=join("|",($c,$s,$e,$o)); for $ex (@f2) { $h{$k}->{$g}->{$ex}=1; }} close(IN); } chomp; $f=$_; ($g,$n,$c,$s,$e,$o,$exons)=split(/\t/,$f); @f2=split(/;/,$exons); $k=join("|",($c,$s,$e,$o)); $n_genes_contained=0; @keys=keys %{$h{$k}}; $n_othergenes=scalar(@keys); for $othergene (@keys) { if($othergene eq $g) { $n_othergenes--; next; } $n3=scalar(keys %{$h{$k}->{$othergene}}); next if($n3 > $n); $n2=0; for $exon (@f2) { if(defined($h{$k}->{$othergene}->{$exon})) { $n2++; } }  if($n2 == $n3) { $n_genes_contained++; } } if($n_othergenes == $n_genes_contained) { print "contained\t$f\n"; }' > ${src}.exonstrings.${src_short}.vs.${target_short}.contained
#cut -f 4-7 ${target}.exonstrings.${src_short}.vs.${target_short}.contained > ${target}.exonstrings.${src_short}.vs.${target_short}.contained.coords
cut -f 2 ${src}.exonstrings.${src_short}.vs.${target_short}.contained | sort -u | sed 's#^#\^#' | sed 's#$#\t#' > ${src}.exonstrings.${src_short}.vs.${target_short}.contained.genes

#fgrep -f ${target}.exonstrings.${src_short}.vs.${target_short}.contained.coords ${src_short}.vs.${target_short}.diff.${src_short} | sort -u | wc -l
#fgrep -v -f ${target}.exonstrings.${src_short}.vs.${target_short}.contained.coords ${src_short}.vs.${target_short}.diff.${src_short} | sort -u > ${src_short}.vs.${target_short}.diff.${src_short}.notcontained
egrep -v -f  ${src}.exonstrings.${src_short}.vs.${target_short}.contained.genes ${src}.exonstrings.${src_short}.vs.${target_short} | cut -f 1-6 > ${src_short}.vs.${target_short}.diff.${src_short}.notcontained.full

done

cat *.notcontained.full | cut -f1 | sed 's#\.'${src_short}'$##' | sort -u | sed 's#^#\t#' | sed 's#$#|#' > ${src_short}.all_notcontained
#fgrep -f ${src_short}.all_notcontained ../gencodeV26.genes.bed | cut -f 4 > ${src_short}.all_notcontained.genes
fgrep -f ${src_short}.all_notcontained ../gencodeV26.genes.bed | cut -d'|' -f 3 | sort | uniq -c > ${src_short}.all_notcontained.gtypes

#get_good_genes_sanity_check.sh
fgrep ".${src_short}" ../disjoint2exons2genes.bed.wstrand.all.diff > ../disjoint2exons2genes.bed.wstrand.all.diff.${src_short}
cat *.notcontained.full | cut -f1 | sort -u | sed 's#^#\t#' | sed 's#$#\t#' > ${src_short}.all_notcontained2
fgrep -v -f ${src_short}.all_notcontained2 ../disjoint2exons2genes.bed.wstrand.all.diff | cut -f 7-11 | head -1 >> ../good_genes_sanity_check2.txt
