#!/usr/bin/env bash
set -exo pipefail
#assume this has the source bucket credentials (not open data) only
#0) starts with manifest of all studies + runs
#in a loop:
#1) checks for finished studies in S3 against manifest are aren't in a) unifier q already and b) aren't done unifying
#2) adds studies from 1) into unifier q
#3) repeat
export LC_ALL=C
#in seconds
SLEEP_TIME=60
#export NUM_EXPECTED_FILES_PER_SAMPLE=36
#export NUM_EXPECTED_FILES_PER_SAMPLE=44
export NUM_EXPECTED_FILES_PER_SAMPLE=40
export PUMP_Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_pump"
export UNIFIER_Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_unify"
export REGION=$(echo "$UNIFIER_Q" | cut -d'.' -f 2)
#export PUMP_S3_OUTPUT="s3://monorail-batch/pump-outputs"
export PUMP_S3_OUTPUT="s3://recount-opendata/recount3expansion/pump"
export UNIFIER_S3_OUTPUT="s3://recount-opendata/recount3expansion/unifier"
export date=$(date +%Y%m%d_%s)


queue_commands () {
    study=$1
    runsF=$2
    for f in `cat $runsF`; do echo "aws sqs send-message --region $REGION --queue-url $PUMP_Q --message-body \"${f}|${study}\""; done > ${runsF}.enqueue
}

    #aws s3 ls --recursive $study_s3 | tr -s " " $'\t' | cut -f 1,3,4 > ${manifestF}.${date}.notdone.${study}.inS3
check_study () {
    #format: idx<tab>study_id<tab>external(run)_id<tab>
    #can use the same manifest file that's used with monitor_finished_pump_studies.sh since this will ignore the header anyway
    study2runF=$1
    #format: tab delimited, long listing that includes file size in bytes (not first field) of all files associated with study's pump run
    #can be either s3 or local filesystem listing, however the path to the file including it's name MUST be the last column!
    #(this is true in both the local POSIX and S3 listings)
    listOfPumpFilesForStudiesF=$2
    #study_id
    study=$3

    return_msg="${study}|DONE"
    #-----First: count number of non-0 byte manifest files as a first check on all samples being finished for the study
    num_samples_expected=$(fgrep $'\t'"$study"$'\t' $study2runF | cut -f 2 | sort -u | wc -l)
    #get full set of pump output files that ARE NOT 0 SIZED for study (all samples run)
    #don't assume we know the # of columns in this file, only assume that the file path is the last
    fgrep "!${study}!" $listOfPumpFilesForStudiesF | tr -s " " $'\t' | fgrep -v $'\t0\t' | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f,-1); $f0=pop(@f); @f2=split(/\//,$f0,-1); $f0=pop(@f2); print "$f0\n";' > ${study}.${date}.finished_files
    #now get just the filenames for those files
    #for fs in `cat ${study}.${date}.finished_files`; do basename $fs; done > ${study}.${date}.finished_files.files
    #now get just the manifest files (1 per sample)
    fgrep ".manifest" ${study}.${date}.finished_files | cut -d'!' -f 1 | sort -u | sed 's#$#\t#' | sed 's#^#\t#' > ${study}.${date}.manifests.samples
    #use that to count the number of samples from the study with non-0 manifest files
    num_samples_finished=$(cat ${study}.${date}.manifests.samples | wc -l)
    #get number of expected runs for this study
    num_samples_expected=$(fgrep $'\t'"$study"$'\t' $study2runF | wc -l)
    if [[ $num_samples_expected -ne $num_samples_finished ]]; then
        #requeue the samples which are missing
        fgrep $'\t'"$study"$'\t' $study2runF | fgrep -v -f ${study}.${date}.manifests.samples | cut -f 3 >${study}.${date}.manifests.samples2rerun
        queue_commands "$study" "${study}.${date}.manifests.samples2rerun"
        return_msg="${study}|MISSING_MANIFESTS"
    else
        #-----Second: do another check on total number of files since infrequently the manifest count check can miss a true failure
        #cut -d'!' -f 1 | sort | uniq -c  
        cut -d'!' -f 1 ${study}.${date}.finished_files | sort | uniq -c | sed 's#^# #' | tr -s " " $'\t' | fgrep $'\t'${NUM_EXPECTED_FILES_PER_SAMPLE}$'\t' | cut -f 3 | sed 's#^#\t#' | sed 's#$#\t#' | sort -u > ${study}.${date}.samples_finished
        
        num_samples_w_expected_num_files=$(cat ${study}.${date}.samples_finished | wc -l)
        if [[ $num_samples_expected -ne $num_samples_w_expected_num_files ]]; then
            #requeue_samples "${manifestF}.${date}.notdone.${study}.finished" "${manifestF}.${date}.notdone.${study}"
            cut -d'!' -f 1 ${study}.${date}.finished_files | sort | uniq -c | sed 's#^# #' | tr -s " " $'\t' | fgrep -v $'\t'${NUM_EXPECTED_FILES_PER_SAMPLE}$'\t' | cut -f 3 | sort -u > ${study}.${date}.samples2rerun
            queue_commands "$study" "${study}.${date}.samples2rerun"
            return_msg="${study}|MISSING_FILES"
        fi
    fi
    echo "$return_msg"
}

study2runF=$1
listOfPumpFilesForStudiesF=$2
study=$3

msg=$(check_study "$study2runF" "$listOfPumpFilesForStudiesF" "$study")
echo "$msg"
