#!/usr/bin/env bash
set -exo pipefail
#0) starts with manifest of all studies + runs
#in a loop:
#1) checks for finished studies in S3 against manifest are aren't in a) unifier q already and b) aren't done unifying
#2) adds studies from 1) into unifier q
#3) repeat
export LC_ALL=C
#export NUM_EXPECTED_FILES_PER_SAMPLE=36
export NUM_EXPECTED_FILES_PER_SAMPLE=44
export PUMP_Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_pump"
export UNIFIER_Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_unify"
export REGION=$(echo "$UNIFIER_Q" | cut -d'.' -f 2)
export PUMP_S3_OUTPUT="s3://monorail-batch/pump-outputs"
export date=$(date +%Y%m%d_%s)

requeue_samples () {
    #e.g. ${manifestF}.${date}.notdone.${study}
    totalF=$1
    #e.g. ${manifestF}.${date}.notdone.${study}.finished
    finishedF=$2
    fgrep -v -f $finishedF $totalF > ${totalF}.samples2requeue
    #requeue at this point or wait for queue to unhide OR push into DLQ?
}

#list of all idx<TAB>studies<TAB>runs being processed by Monorail Pump
manifestF=$1
#optional, list of studies which are finished pump processing, <TAB>studies<TAB>
doneF=$2

if [[ -z $doneF ]]; then
    doneF="STUDIES_QUEUED_FOR_UNIFIER.txt" 
fi 

wc=0
if [[ -e $doneF ]]; then
    wc=$(cat $doneF | wc -l)
fi
if [[ $wc -gt 0 ]]; then
    fgrep -v -f $doneF $manifestF > ${manifestF}.${date}.notdone
else
    ln -fs $manifestF ${manifestF}.${date}.notdone
fi

cut -f 2 ${manifestF}.${date}.notdone | sed 's#^#/#' | sed 's#$#/#' | sort -u > ${manifestF}.${date}.notdone.studies

echo -n "" > s3.pump.${date}
top_levels=$(aws s3 ls $PUMP_S3_OUTPUT/ | tr -s " " $'\t' | cut -f3)
for d in $top_levels; do
    aws s3 ls $PUMP_S3_OUTPUT/$d | tr -s " " $'\t' | cut -f 3 | sed 's#^#'$PUMP_S3_OUTPUT'/'$d'#' >> s3.pump.${date}
done
fgrep -f ${manifestF}.${date}.notdone.studies s3.pump.${date} > s3.pump.${date}.notdone

wc=$(cat s3.pump.${date}.notdone | wc -l)
if [[ $wc -gt 0 ]]; then
   for study_s3 in `cat s3.pump.${date}.notdone`; do
        study=$(basename $study_s3 | sed 's#/$##')
        #-----First: count number of non-0 byte manifest files as a first check on all samples being finished for the study
        fgrep $'\t'"$study"$'\t' ${manifestF}.${date}.notdone | sed 's#$#\t#' > ${manifestF}.${date}.notdone.${study}
        num_samples_expected=$(cat ${manifestF}.${date}.notdone.${study} | wc -l)
        aws s3 ls --recursive $study_s3 | tr -s " " $'\t' | cut -f 1,3,4 > ${manifestF}.${date}.notdone.${study}.inS3
        cat ${manifestF}.${date}.notdone.${study}.inS3 | fgrep manifest | fgrep -v $'\t0\t' > ${manifestF}.${date}.notdone.${study}.manifestsINs3
        num_samples_finished=$(cut -f 3 ${manifestF}.${date}.notdone.${study}.manifestsINs3 | cut -d'.' -f 1 | sort -u | wc -l)
        if [[ $num_samples_expected -ne $num_samples_finished ]]; then
            #requeue the samples which are missing
            echo -n "" > ${manifestF}.${date}.notdone.${study}.manifestsINs3.finished
            for sample in `cut -f 3 ${manifestF}.${date}.notdone.${study}.manifestsINs3 | cut -d'.' -f 1 | sort -u`; do 
                basename $sample | sed 's#^#\t#' | sed 's#$#\t#' >> ${manifestF}.${date}.notdone.${study}.finished
            done
            requeue_samples "${manifestF}.${date}.notdone.${study}.finished" "${manifestF}.${date}.notdone.${study}"
        else
            #-----Second: do another check on total number of files since infrequently the manifest count check can miss a true failure
            cat ${manifestF}.${date}.notdone.${study}.inS3 | cut -f 3 | cut -d'/' -f 1-5 | sort | uniq -c | tr -s " " $'\t' | fgrep $'\t'${NUM_EXPECTED_FILES_PER_SAMPLE}$'\t' | cut -d'.' -f 1 | cut -f 3 | cut -d'/' -f 5 | sed 's#^#\t#' | sed 's#$#\t#' | sort -u > ${manifestF}.${date}.notdone.${study}.finished
            num_samples_w_expected_num_files=$(cat ${manifestF}.${date}.notdone.${study}.finished | wc -l)
            if [[ $num_samples_expected -ne $num_samples_w_expected_num_files ]]; then
                #requeue_samples "${manifestF}.${date}.notdone.${study}.finished" "${manifestF}.${date}.notdone.${study}"
                echo "NOT_DONE	$study	skipping"
            else
                #----Third: all samples done, queue up for Unifier
                study_s3=$(echo "$study_s3" | sed 's#/$##')
                aws sqs send-message --region $REGION --queue-url $UNIFIER_Q --message-body "$study_s3"
                echo $'\t'"$study"$'\t' >> $doneF
            fi
        fi
    done
fi
