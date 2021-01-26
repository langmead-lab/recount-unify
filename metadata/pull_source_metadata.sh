#!/usr/bin/env bash
#use this to get SRA sourced, per-sample metadata (potentially others to be supported in future)

#e.g. ERP001942
study=$1
#hg38 or grcm38
ref=$2
#0 (no) or 1 (yes) to include dbGaP in the query
protected=$3
#optional, if set we only for query for runs from the the study which are labeled as "transcriptomic"
#skips if defined
skip_nontranscriptomic=$4

orgn="human"
if [[ "$ref" == "grcm38" ]]; then
   orgn="mouse"
fi

dir=$(dirname $0)
extra=
if [[ -n $protected && $protected -gt 0 ]]; then
    extra="--include-protected"
fi

#replace spaces with "_" to allow for use in filenames
orgn_orig=$orgn
orgn=`echo -n "$orgn_orig" | perl -ne '$o=$_; $o=~s/\s+/_/g; print "$o";'`
#assumes write to the current working directory
mkdir -p xml_out
mkdir -p err_out

python $dir/../recount-pump/metadata/scripts/fetch_sra_metadata.py --accession $study --orgn $orgn --xml-path xml_out --err-path err_out $extra > fetch_sra_${orgn}.txt 2>&1

/bin/bash -x fetch_${orgn}.jobs > fetch_${orgn}.jobs.run 2>&1
/bin/bash -x parse_${orgn}.sh > parse_${orgn}.sh.run 2>&1

#try to find non-transcriptomic runs from this study, if not skipped 
if [[ -z $skip_nontranscriptomic ]]; then
    #pick up any runs which don't have "transcriptomic" as source but share a study with ones that do
    mkdir nontranscriptomic
    pushd nontranscriptomic
    python $dir/../recount-pump/metadata/scripts/fetch_sra_metadata.py --accession $study --orgn $orgn --xml-path xml_out --err-path err_out $extra --non-transcriptomic > fetch_sra_${orgn}.txt 2>&1
    no_records=$(fgrep -m1 'Total # of records is 0 ' fetch_sra_${orgn}.txt)

    if [[ -z $no_records ]]; then
        /bin/bash -x fetch_${orgn}.jobs > fetch_${orgn}.jobs.run 2>&1
        /bin/bash -x parse_${orgn}.sh > parse_${orgn}.sh.run 2>&1
        cut -f 2 all_${orgn}_sra.tsv | sort -u | sed -e 's/$/\t/' > all_${orgn}_sra.tsv.studies
        fgrep -f all_${orgn}_sra.tsv.studies ../all_${orgn}_sra.tsv | cut -f 2 | sed -e 's/$/\t/' | sort -u > studies_found_in_transcriptomic
        fgrep -f studies_found_in_transcriptomic all_${orgn}_sra.tsv > all_${orgn}_sra.tsv.in_transcriptomic
        cat ../all_${orgn}_sra.tsv all_${orgn}_sra.tsv.in_transcriptomic | sort -u > ../all_${orgn}_sra.with_nontranscriptomic_runs.tsv
        popd
        ln -fs all_${orgn}_sra.with_nontranscriptomic_runs.tsv all.runs.tsv
    else
        popd
        ln -fs all_${orgn}_sra.tsv all.runs.tsv
    fi
else
    ln -fs all_${orgn}_sra.tsv all.runs.tsv
fi

#get study level info
cut -f 2,8-10 all.runs.tsv | sort -u > all.runs.tsv.study_level

#add header
cat <(echo "run_acc	study_acc	sample_acc	experiment_acc	submission_acc	submission_center	submission_lab	study_title	study_abstract	study_description	experiment_title	design_description	sample_description	library_name	library_strategy	library_source	library_selection	library_layout	paired_nominal_length	paired_nominal_stdev	library_construction_protocol	platform_model	sample_attributes	experiment_attributes	spot_length	taxon_id	sci_name	common_name	sample_name	sample_title	sample_bases	sample_spots	run_published	size	run_total_bases	run_total_spots	num_reads	num_spots	read_info	run_alias	run_center_name	run_broker_name	run_center	inferred_read_length	inferred_total_read_count") all.runs.tsv > all.runs.tsv.1
mv all.runs.tsv.1 all.runs.tsv
