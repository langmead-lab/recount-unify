#!/usr/bin/env bash
#use this to get SRA sourced, per-sample metadata (potentially others to be supported in future)

#e.g. ERP001942
study=$1
#hg38 or grcm38
ref=$2
#the recount project level metadata for the study, so we can map the rail_ids
#e.g. sra.recount_project.ERP001942.MD.gz
rc_project_file=$3
#e.g. sra
dsource=$4

perl -e '$study="'$study'"; $study=~/(..)$/; $lo=$1; `mkdir -p metadata/$lo/$study`; print "metadata/$lo/$study\n";' > ${study}.dir
sdir=`cat ${study}.dir`

orgn="human"
if [[ "$ref" == "grcm38" ]]; then
   orgn="mouse"
fi

dir=$(dirname $0)

#replace spaces with "_" to allow for use in filenames
orgn_orig=$orgn
orgn=`echo -n "$orgn_orig" | perl -ne '$o=$_; $o=~s/\s+/_/g; print "$o";'`

#change into study specific working directory
working_dir=${study}_metadata_temp
mkdir -p $working_dir
pushd $working_dir

#assumes write to the current working directory
mkdir -p xml_out
mkdir -p err_out

if [[ -z $PYTHON_PATH ]]; then
    PYTHON_PATH=python3
fi

#want as relaxed a query pattern as possible, we'll filter in only those runs we are actually processing, the rest will be discarded
$PYTHON_PATH $dir/../recount-pump/metadata/scripts/fetch_sra_metadata.py --base-query "${study}[Accession]" --orgn $orgn --xml-path xml_out --err-path err_out > fetch_sra_${orgn}.txt 2>&1

/bin/bash -x fetch_${orgn}.jobs > fetch_${orgn}.jobs.run 2>&1
/bin/bash -x parse_${orgn}.sh > parse_${orgn}.sh.run 2>&1

ln -fs all_${orgn}_sra.tsv all.runs.tsv

#get study level info
cut -f 2,8-10 all.runs.tsv | sort -u > all.runs.tsv.study_level

#add header
cat <(echo "external_id	study	sample_acc	experiment_acc	submission_acc	submission_center	submission_lab	study_title	study_abstract	study_description	experiment_title	design_description	sample_description	library_name	library_strategy	library_source	library_selection	library_layout	paired_nominal_length	paired_nominal_stdev	library_construction_protocol	platform_model	sample_attributes	experiment_attributes	spot_length	taxon_id	sci_name	common_name	sample_name	sample_title	sample_bases	sample_spots	run_published	size	run_total_bases	run_total_spots	num_reads	num_spots	read_info	run_alias	run_center_name	run_broker_name	run_center	inferred_read_length	inferred_total_read_count") all.runs.tsv > all.runs.tsv.1

#remove taxon_id, sci_name, and common_name as these were added later (post-recount3), as well as the dubious and potentially redundant inferred_read_length & inferred_total_read_count columns at the end
cat all.runs.tsv.1 | cut -f 1-25,29-43 > all.runs.tsv

#now map in rail_ids
cat <(gzip -cd ../$rc_project_file | cut -f 1-3) all.runs.tsv | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f,-1); if(scalar(@f) == 3) { $hmap{$f[1]."\t".$f[2]}=$f[0]; next; } $k=$f[0]."\t".$f[1]; $rail_id=$hmap{$k}; if(!defined($rail_id)) { print STDERR "couldnt map $k to a rail_id, skipping\n"; next; } print "$rail_id\t$f\n";' > all.runs.tsv.2
cat all.runs.tsv.2 | gzip > ../$sdir/${dsource}.${dsource}.${study}.MD.gz

#cleanup
popd
#rm -rf $working_dir ${study}.dir
#rm -f all.runs.tsv.1 all_${orgn}_sra.tsv all_${orgn}_sra.with_nontranscriptomic_runs.tsv all.runs.tsv.? ../${study}.dir
