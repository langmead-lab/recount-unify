date="2019-10-01"
data_source="data_sources"
organism="Homo sapiens"

sdir=$(dirname $0)

#assumes we're using samples.tsv and samples.tsv.header
study=$1
#human, mouse, gtex, or tcga
org=$2
#sra, gtex, or tcga
dsource=$3
#RNA-seq, tissue, celltype Predictions/curations file
#e.g. Human.SRA_curated_predict.meta.w_rids.tsv
#or
#Mouse.SRA_curated_predict.meta.w_rids.w_cell_types_nas.tsv
#not used with gtex or tcga
preds_file=$4

perl -e '$study="'$study'"; $study=~/(..)$/; $lo=$1; `mkdir -p $lo/$study`; print "$lo/$study\n";' > ${study}.dir
dir=`cat ${study}.dir`

#head -1 samples.tsv > samples.tsv.header
cat $sdir/samples.tsv.header <(fgrep "	$study	" $sdir/samples.tsv) > $dir/samples.tsv.${study}
cut_cols='1-3'
if [[ "$org" == "mouse" ]]; then
    organism="Mus musculus"
    date="2020-01-01"
    #sra specific MD file
    cut -f 1-43 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.${dsource}.${study}.MD.gz
    #QC MD file
    cut -f 1-3,44-152 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.recount_qc.${study}.MD.gz
fi
if [[ "$org" == "human" ]]; then
    #sra specific MD file
    cut -f 1-41 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.${dsource}.${study}.MD.gz
    #QC MD file
    cut -f 1-3,42-150 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.recount_qc.${study}.MD.gz
fi
if [[ "$org" == "gtex" ]]; then
    date="2019-11-01"
    #sra specific MD file
    cut -f 1-71 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.${dsource}.${study}.MD.gz
    #QC MD file
    cut -f 1,3,4,72-180 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.recount_qc.${study}.MD.gz
    cut_cols='1,3,4'
fi
if [[ "$org" == "tcga" ]]; then
    date="2019-12-01"
    #sra specific MD file
    cut -f 1-840 $dir/samples.tsv.${study} | gzip > $dir/${dsource}.${dsource}.${study}.MD.gz
    #QC MD file
    cut -f 1-3,841- $dir/samples.tsv.${study} | gzip > $dir/${dsource}.recount_qc.${study}.MD.gz
fi


proj_header='rail_id	external_id	study	project	organism	file_source	metadata_source	date_processed'
dsource_path="$data_source/$dsource"
extra_proj_cols="$organism	$dsource_path	$dsource_path	$date"

#project MD file
cat <(echo "$proj_header") <(tail -n+2 $dir/samples.tsv.${study} | cut -f $cut_cols | perl -ne 'chomp; $f=$_; ($rid,$run,$study)=split(/\t/,$f); print "$f\t$study\t'"$extra_proj_cols"'\n";') | gzip > $dir/${dsource}.recount_project.${study}.MD.gz
#pcat $dir/${dsource}.recount_project.${study}.MD.gz | tail -n+2 >> all.recout_project.MD


#RNA-seq, tissue, celltype Predictions/curations
#if [[ "$org" == "mouse" || "$org" == "human" ]]; then
if [[ -n $preds_file ]]; then
    cat <(head -1 $preds_file) <(fgrep "	$study	" $preds_file) | gzip > $dir/${dsource}.recount_pred.${study}.MD.gz
fi

#add seqtk in, assume there's a generically named file with the normal join key triplet (rail_id, external_id, and study) prepended + the actual seqTK QC fields file named "seqtk.w_study.tsv.cut"
lo=$(perl -e '$s="'$study'"; $s=~/(..)$/; $lo=$1; print "$lo\n";')
cat <(head -1 $sdir/seqtk.w_study.tsv.cut) <(fgrep "	$study	" $sdir/seqtk.w_study.tsv.cut) | gzip > ../metadata/${lo}/${study}/sra.recount_seq_qc.${study}.MD.gz

rm $dir/samples.tsv.${study}
rm ${study}.dir
