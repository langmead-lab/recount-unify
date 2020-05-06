#!/usr/bin/env bash
#once all tranches of a data_source are finished
#this script is run to to collect (and split) by A) study and B) annotation
#the 1) gene sums and 2) exon sums

dir=$(dirname $0)

#e.g. "srav3_human" or "srav1_mouse"
tranche_prefix=$1
#e.g. for SRAv3Human and SRAv1Mouse this is 9
max_tranche_num=$2
#e.g. G029.G026.R109.F006.20190220.gtf.exons2genes
exonids2annotations=$3
#e.g. G026,G029,R109,F006,ERCC,SIRV
annotations=$4
#e.g. all.exon_counts.rejoined.tsv.gz.coords, generated by the per tranche snakemake, but only need one (any tranche will work)
rejoined_exon_coords=$5

#number of parallel processors to use, on the elephant machines: ~40
num_procs=$6

#e.g. G026
main_annotation=`echo "$annotations" | cut -d',' -f 1`

#e.g. 1709834 for SRAv3Human
num_exon_rows=`cat $rejoined_exon_coords | wc -l`

#get fullpath to output dir so we can write to it from anywhere
outdir=`pwd`/${tranche_prefix}.exon_splits_per_study
mkdir -p $outdir
        
#run exons splits into per-study files (split jobs file is already created in the per-tranche Snakefile run)
for i in {0..${max_tranche_num}}; do
    mkdir -p ${outdir}/tranche${i}
    tdir=${tranche_prefix}${max_tranche_num}
    pushd $tdir
    pypy ${dir}/create_exon_sums_by_study_splits.py ids.tsv all.exon_counts.rejoined.tsv.gz all.exon_counts.rejoined.tsv.gz.accession_header ${outdir}/tranche${i} > exon_sums.splits.${i}.jobs
    /usr/bin/time -v parallel -j ${num_procs} < exon_sums.splits.${i}.jobs > exon_sums.splits.${i}.jobs.run 2>&1
    popd
done
    
bitmasks_file=exon_counts.bitmasks
 
#create bitmasks so we can quickly split exon sums already split into per-study files further into per-annotation files
cat ${rejoined_exon_coords} | pypy ${dir}/../rejoin/map_exon_sum_rows_to_annotations.py $exonids2annotations $annotations > ${bitmasks_file}.tsv

bitmask_coords_file=${bitmasks_file}.coords
#make 1-based, not BED
cut -f 3-5 ${bitmasks_file}.tsv | perl -ne 'chomp; ($c,$s,$e)=split(/\t/,$_); $s++; print "$c|$s|$e\n";' > ${bitmask_coords_file}

#create per-annotation split jobs;
echo -n "" > exon_sums.per_annotation.splits.jobs
for sfile in `ls ${outdir}/tranche?/*.*.gz`; do
    study=`basename $sfile | cut -d'.' -f 2`
    lo=`echo $study | sed 's/^.*\(..\)$/\1/'`
    #for purely organizational reasons, we make subdirs from the 2 lowest-ordER characters of the study accession
    mkdir -p exon_splits_done/$lo
    final_outdir=exon_splits_done/$lo
    #/bin/bash -x split_out_exons.sh G026,G029,R109,ERCC,SIRV all.exon_counts.rejoined.tsv.gz.coords.bitmasks2 1709834 ERP001942 examples_done/ERP001942.gz examples_done/ all.exon_counts.rejoined.tsv.gz.coords.bitmasks2_1base_coords > ERP001942.final_annotation_split.run 2>&1
    echo "/usr/bin/time -v /bin/bash -x ${dir}/../rejoin/split_out_exons.sh \"${annotations}\" ${bitmasks_file}.tsv $num_exon_rows $study $sfile $final_outdir ${bitmask_coords_file} > $final_outdir/${study}.annot_split.run 2>&1" >> exon_sums.per_annotation.splits.jobs
done 

/usr/bin/time -v parallel -j ${num_procs} < exon_sums.per_annotation.splits.jobs > exon_sums.per_annotation.splits.jobs.run 2>&1
