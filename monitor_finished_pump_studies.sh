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
export UNIFIER_Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_unify"
export REGION=$(echo "$UNIFIER_Q" | cut -d'.' -f 2)
#export PUMP_S3_OUTPUT="s3://monorail-batch/pump-outputs"
export PUMP_S3_OUTPUT="s3://recount-opendata/recount3expansion/pump"
export UNIFIER_S3_OUTPUT="s3://recount-opendata/recount3expansion/unifier"
export date=$(date +%Y%m%d_%s)

#list of all idx<TAB>studies<TAB>runs being processed by Monorail Pump
manifestF=$1
#human or mouse
org0=$2
if [[ -z $org0 ]]; then
    org0="human"
fi

#for fgrepping later
head -1 $manifestF | tr $'\t' $'\n' | fgrep -n "" > ${manifestF}.fields
study_col=$(egrep -e $':study_id$' ${manifestF}.fields | cut -d':' -f 1)
sample_col=$(egrep -e $':external_id$' ${manifestF}.fields | cut -d':' -f 1)
if [[ $sample_col -lt $study_col ]]; then
    paste <(cut -f $study_col $manifestF) <(cut -f $sample_col $manifestF) | sed 's#^#\t' | sed 's#$#\t#' > ${manifestF}.cut.tabs
else
    cut -f ${study_col},${sample_col} $manifestF | sed 's#^#\t' | sed 's#$#\t#' > ${manifestF}.cut.tabs
fi
mkdir -p PUMP_DONES
export PUMP_DONES=$PUMP_S3_OUTPUT/$org0/DONES

while true; do
    aws s3 sync --no-sign-request $PUMP_DONES/ ./PUMP_DONES/
    #example of path in the .DONE file:
    #ERP001942.ERR188246.20240807_1723047815.DONE
    #containing:
    #s3://recount-opendata/recount3expansion/pump/human/95/SRP222095/69/SRR10133969.20240807_1723052789

    #determine which studies are finished (ignore unfinished studies for now, another script will take care of stragglers):
    pushd ./PUMP_DONES/
    #get study,run ids for the DONE samples
    ls *.DONE | cut -d'.' -f 1-2 | sort -u > dones
    #count number of finished run ids for each study
    cut -f 1 dones | sort | uniq -c | sed 's#^# #' | tr -s " " $'\t' | cut -f 2- | sort > dones.counts2study
    #format the list of studies for grepping
    cut -f 2 dones.counts2study | sed 's#$#\t#' | sed 's#^#\t#' > dones.counts2study.tabs 
    #grep out the list of studies and get the expected run counts for each study to compare with
    fgrep -f dones.counts2study.tabs ${manifestF}.cut.tabs | cut -f 2-3 | sort | uniq -c | sed 's#^# #' | tr -s " " $'\t' | cut -f 2- | sort > ${manifestF}.cut.tabs.subset_counts

    #compare finished # of runs per study counts with expected counts, get the studies that match
    comm -1 -2 ${manifestF}.cut.tabs.counts2study dones.counts2study > dones.counts2study.matching
    #get just the study ids of the matching ones
    cut -f 2 dones.counts2study.matching > dones.counts2study.matching.studies

    #now queue up finished studies from pump, for unifying
    for study0 in `cat dones.counts2study.matching.studies`; do
        lo=${study0: -2}
        #send this to the queue
        studys3="s3://recount-opendata/recount3expansion/pump/$org0/$lo/$study0"
        echo "$studys3" >> ../ready2unify.studies 
        #clear pump .DONE files from S3 for this study
        #the idea being to try to avoid/minimize as best as possible any multiple overlapping enqueues of the same study
        #(only ever want 1!)
        for f in `ls ${study0}.*`; do
            aws s3 rm $PUMP_S3_OUTPUT/$org0/DONES/$f
            #still track these pump files just in case...
            aws s3 cp $f $PUMP_S3_OUTPUT/$org0/UNIFYING/
        done
        #clear local directory of sample .DONE files for this study as well
        rm -rf ${study0}.*
        #finally, when all other sentinel files have been (re)moved, enqueue the unifier
        aws sqs send-message --region $REGION --queue-url $UNIFIER_Q --message-body "$studys3"
    done
    popd
    sleep $SLEEP_TIME
done

#mkdir -p UNIFIER_DONES
#export UNIFIER_DONES=$UNIFIER_S3_OUTPUT/$org0/DONES
#aws s3 sync --no-sign-request $UNIFIER_DONES/ ./UNIFIER_DONES/
