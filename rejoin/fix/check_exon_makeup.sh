#e.g. G026
src=$1
#e.g. G029
target=$2

export LC_ALL=C

fgrep "$src" ../disjoint2exons2genes.bed > disjoint2exons2genes.bed.${src}
fgrep "$target" ../disjoint2exons2genes.bed > disjoint2exons2genes.bed.${target}

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

diff ${src}.exonstrings.bycoords ${target}.exonstrings.by${src_short} > ${src_short}.vs.${target_short}.diff
fgrep "<" ${src_short}.vs.${target_short}.diff | cut -d' ' -f 2- | cut -f 1-4 > ${src_short}.vs.${target_short}.diff.${src_short}

fgrep -f ${src_short}.vs.${target_short}.diff.${src_short} ${src}.exonstrings | sort -u | sort -k3,3 -k4,4n -k5,5n > ${src}.exonstrings.${src_short}.vs.${target_short}
fgrep -f ${src_short}.vs.${target_short}.diff.${src_short} ${target}.exonstrings | sort -u | sort -k3,3 -k4,4n -k5,5n > ${target}.exonstrings.${src_short}.vs.${target_short}

cat ${target}.exonstrings.${src_short}.vs.${target_short} | perl -ne 'BEGIN { open(IN,"<'${src}.exonstrings.${src_short}.vs.${target_short}'"); while($line=<IN>) { chomp($line); ($g,$n,$c,$s,$e,$o,$exons)=split(/\t/,$line); @f2=split(/;/,$exons); $k=join("|",($c,$s,$e,$o)); for $ex (@f2) { $h{$k}->{$ex}=1; }} close(IN); } chomp; $f=$_; ($g,$n,$c,$s,$e,$o,$exons)=split(/\t/,$f); @f2=split(/;/,$exons); $k=join("|",($c,$s,$e,$o)); $n2=0; for $ex (@f2) { if(defined($h{$k}->{$ex})) { $n2++; } } if($n2 == $n) { print "contained\t$f\n"; }' > ${target}.exonstrings.${src_short}.vs.${target_short}.contained
cut -f 4-7 ${target}.exonstrings.${src_short}.vs.${target_short}.contained > ${target}.exonstrings.${src_short}.vs.${target_short}.contained.coords

#fgrep -f ${target}.exonstrings.${src_short}.vs.${target_short}.contained.coords ${src_short}.vs.${target_short}.diff.${src_short} | sort -u | wc -l
fgrep -v -f ${target}.exonstrings.${src_short}.vs.${target_short}.contained.coords ${src_short}.vs.${target_short}.diff.${src_short} | sort -u > ${src_short}.vs.${target_short}.diff.${src_short}.notcontained

cat *.notcontained | sort -u | cut -f 1-3 | sed 's#$#\t#' > ${src_short}.all_notcontained
fgrep -f ${src_short}.all_notcontained ../gencodeV26.genes.bed | cut -f 4 > ${src_short}.all_notcontained.genes
cut -d'|' -f 3 ${src_short}.all_notcontained.genes | sort | uniq -c > ${src_short}.all_notcontained.genes.types
