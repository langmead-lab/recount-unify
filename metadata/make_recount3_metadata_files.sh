#!/usr/bin/env bash
#create 2 of the 3 required metadata files required for a recount3 loadable study
#(the 3rd is the source metadata, which is handled by pull_source_metadata.sh)
#1) recount project (manifest of samples and their data sources/organisms mainly for loading)
#2) original sample metadata (either from SRA or custom generated by the source organization)

#e.g. 2021-01-21
date=$(date +%Y-%m-%d)
data_source="data_sources"
organism="Homo sapiens"

sdir=$(dirname $0)

study=$1
#hg38 or grcm38
org=$2
#assumes project qc metadata has already been finalized (but NOT gzipped yet)
qc_file=$3
#"sra" or custom datasource
dsource=$4

if [[ -z $dsource ]]; then
    dsource="sra"
fi

perl -e '$study="'$study'"; $study=~/(..)$/; $lo=$1; `mkdir -p metadata/$lo/$study`; print "metadata/$lo/$study\n";' > ${study}.dir
dir=`cat ${study}.dir`

#1) recount project format
#rail_id external_id     study   project organism        file_source     metadata_source date_processed
#313542  SRR10230231     SRP224299       SRP224299       Homo sapiens    data_sources/sra        data_sources/sra        2019-10-01
if [[ "$org" == "grcm38" ]]; then
    organism="Mus musculus"
fi

proj_header='rail_id	external_id	study	project	organism	file_source	metadata_source	date_processed'
dsource_path="$data_source/$dsource"

#create recount_project metadata file
cat <(echo "$proj_header") <(cat $qc_file | tail -n+2 | fgrep "$study	" | cut -f 1-3 | perl -ne 'chomp; $f=$_; print "$f\t'$study'\t'"$organism"'\t'$dsource_path'\t'$dsource_path'\t'$date'\n";') | gzip > $dir/${dsource}.recount_project.${study}.MD.gz

echo "RC=$dir/${dsource}.recount_project.${study}.MD.gz"

#2) original sample metadata can be *any* tab delimited set of fields as long as it starts with these 3 columns for each row:
#rail_id external_id     study
num_cols=$(head -1 samples.tsv | sed 's/rail_id\trun\t/rail_id\texternal_id\t/' | tee samples.tsv.header | tr \\t \\n | wc -l)
num_cols_minus_jxs=$(( num_cols - 3 ))
cat <(cut -f 1-${num_cols_minus_jxs} samples.tsv.header) <(fgrep "$study	" samples.tsv | cut -f 1-${num_cols_minus_jxs}) | gzip > $dir/${dsource}.${dsource}.${study}.MD.gz

echo "RPROJ=$dir/${dsource}.${dsource}.${study}.MD.gz"

#3) just grep out this study's samples & gzip the finalized QC file into place
cat <(head -1 $qc_file) <(fgrep "$study	" $qc_file) | gzip > $dir/${dsource}.recount_qc.${study}.MD.gz

echo "QC=$dir/${dsource}.recount_qc.${study}.MD.gz"

rm ${study}.dir
