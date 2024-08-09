#!/usr/bin/env bash
set -exo pipefail
#assume the AWS Open Data Creds are available and are being used for these operations
export LOG_DIR=../transfer_logs_unifier
mkdir -p $LOG_DIR
#in seconds
SLEEP_TIME=60
#source bucket
SB="s3://monorail-batch"
#target bucket
TB="s3://recount-opendata"
awsprofile=$1
if [[ -z $awsprofile ]]; then
    awsprofile="default"
fi
#human or mouse
org0=$2
if [[ -z $org0 ]]; then
    org0="human"
fi
pSB=$SB/unifier-outputs/$org0
pTB=$TB/recount3expansion/unifier/$org0
pDONES="s3://monorail-batch/UNIFIER_DONES"

while true; do
    #remove any previous done files to avoid confusion about what's already been copied
    rm -rf *.DONE
    #aws s3 ls $pDONES | fgrep ".DONE" | tr -s " " $'\t' | cut -f 4 > pdones
    #get latest set of finsihed samples (from pump)
    aws s3 --profile $awsprofile cp --request-payer requester --recursive $pDONES/ ./
    #for each sample (.DONE file):
    for f in `ls *.DONE`; do
        #get the full S3 source path to copy the pump results from
        #e.g. s3://monorail-batch/unifier-output/human/metadata/95/SRP222095
        s3path=$(cat $f)
        suffix=$(echo "$s3path" | perl -ne 'chomp; $f=$_; @f=split(/\//,$f,-1); $study=pop(@f); $lo=pop(@f); print join("/",($lo,$study))."\n";')
        #use this for local logging purposes
        logd=$LOG_DIR/$suffix
        mkdir -p $logd
        #copy back the full pump results directory to the target bucket
        /usr/bin/time -v aws s3 --profile $awsprofile sync --request-payer requester $s3path/ $pTB/$suffix/ > $logd/s3copy 2>&1
        /usr/bin/time -v aws s3 --profile $awsprofile cp --request-payer requester $pDONES/$f $pTB/DONES/ > $logd/s3copy.done 2>&1
        #then remove path from source bucket
        /usr/bin/time -v aws s3 --profile $awsprofile rm --request-payer requester --recursive $s3path/ > $logd/s3rm 2>&1
        /usr/bin/time -v aws s3 --profile $awsprofile rm --request-payer requester $pDONES/$f > $logd/s3rm.done 2>&1
    done
    sleep $SLEEP_TIME
done
