#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit 
#assumes zstd, tabix, bgzip, and GNU parallel are in PATH

d=$(dirname $0)

#GENOME_FA='hg38_plus.fa'
#GENOME_FA_SIZES='hg38_plus.fa.header.sizes.first25'

GENOME_FA='mm10.fa'
GENOME_FA_SIZES='mm10.fa.header.sizes'

#find -L ./ERP015294 -name "*.jx_bed.zst" > ERP015294.all_jx_bed_files
#egrep -e '_att0' ERP015294.all_jx_bed_files > ERP015294.all_jx_bed_files.done
#egrep -v -e 'umapped' ERP015294.all_jx_bed_files.done > ERP015294.all_jx_bed_files.done.mapped

#path to file with list of paths to zstd compressed regtools-extract jx files
mapped_jx_zstds=$1
study=$2
num_procs=$3
prefix=$study

rm -rf ${study}_regtools
mkdir ${study}_regtools
#have to work around the fact that perl tries to interpolate the "@" symbols in the path
extract="${d}/extract_and_consolidate_regtools_jxs_per_sample.sh"
#ln -fs $d/${extract} ${study}_regtools/
#cat $mapped_jx_zstds | perl -ne 'BEGIN { $p="'${study}_regtools'"; } chomp; $f=$_; $f1=$f; @f=split(/\//,$f); $f1=pop(@f); $f1=~s/\.zst$//; $f=~s/!/\\!/g; $f1=~s/!/_/g; print "zstd -cd $f | $p/extract_and_consolidate_regtools_jxs_per_sample.sh > $p/$f1.condensed 2>$p/$f1.err\n";' > ${prefix}.zstd_jobs
cat $mapped_jx_zstds | perl -ne 'BEGIN { $p="'${study}_regtools'"; } chomp; $f=$_; $f1=$f; @f=split(/\//,$f); $f1=pop(@f); $f1=~s/\.zst$//; $f=~s/!/\\!/g; $f1=~s/!/_/g; print "zstd -cd $f | /bin/bash '$extract' > $p/$f1.condensed 2>$p/$f1.err\n";' > ${prefix}.zstd_jobs

parallel -j $num_procs < ${prefix}.zstd_jobs

#get actual jx coords; collapse ones with same coords (after sorting by coords); cut to proper input format for merge.py
ls ${study}_regtools/*.condensed | perl -ne 'BEGIN { $i=0; } chomp; print "$_\t".($i++)."\n";' > ${prefix}.condensed.merge.manifest

#merge jx's from individual samples into one set tracking read counts per sample
time pypy ${d}/../merge/merge.py --list-file ${prefix}.condensed.merge.manifest > ${prefix}.condensed.merge.manifest.merged 2>${prefix}.condensed.merge.manifest.err

#extract actual motifs for the jxs
time cat ${prefix}.condensed.merge.manifest.merged | ${d}/perbase -c $GENOME_FA_SIZES -g $GENOME_FA -f $GENOME_FA_SIZES > ${prefix}.condensed.merge.manifest.merged.tsv.motifs 2>errs

paste <(cut -f 1-4 ${prefix}.condensed.merge.manifest.merged.tsv.motifs) <(cut -f 7 ${prefix}.condensed.merge.manifest.merged.tsv.motifs) <(cut -f 6 ${prefix}.condensed.merge.manifest.merged.tsv.motifs) > ${prefix}.condensed.merge.manifest.merged.tsv.motifs.pasted

#finally annotate the jxs
#time cat ${prefix}.condensed.merge.manifest.merged.tsv.motifs.pasted | pypy ${d}/../annotate/annotate_sjs.py --compiled-annotations ${d}/../annotate/annotation/annotated_junctions.tsv.gz > ${prefix}.condensed.merge.manifest.merged.annotated.tsv

#cat ${prefix}.condensed.merge.manifest.merged.annotated.tsv | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=split(/,/,$f[11]); shift(@f1); $c=scalar(@f1); $count=0; map { ($sid,$r)=split(/:/,$_); $count+=$r; } @f1; print "$f\t$c\t$count\n";' | bgzip > ${prefix}.snaptron_lite.tsv.bgz

#tabix -s 2 -b 3 -e 4 ${prefix}.snaptron_lite.bgz
