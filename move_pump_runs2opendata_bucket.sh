#!/usr/bin/env bash
set -exo pipefail
#assume the AWS Open Data Creds are available and are being used for these operations
export LC_ALL=C
export LOG_DIR=../transfer_logs
mkdir -p $LOG_DIR
export threads=10
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
pSB=$SB/pump-outputs/$org0
pTB=$TB/recount3expansion/pump/$org0
pDONES="s3://monorail-batch/PUMP_DONES"

#REFACTORED to do parallel copyies/removals

while true; do
    #remove any previous done files to avoid confusion about what's already been copied
    rm -rf *.DONE
    #get latest set of finsihed samples (from pump)
    aws s3 --profile $awsprofile cp --request-payer requester --recursive $pDONES/ ./
    echo -n "" > $LOG_DIR/all_copies.${date}
    echo -n "" > copy_finished_pump_runs2opendata.jobs
    #for each sample (.DONE file):
    for f in `ls *.DONE`; do
        #get the full S3 source path to copy the pump results from
        s3path=$(cat $f)
        suffix=$(echo "$s3path" | perl -ne 'chomp; $f=$_; @f=split(/\//,$f,-1); $sdate=pop(@f); $lo2=pop(@f); $study=pop(@f); $lo=pop(@f); print join("/",($lo,$study,$lo2,$sdate))."\n";')
        #use this for local logging purposes
        logd=$LOG_DIR/$suffix
        mkdir -p $logd
        #copy back the full pump results directory to the target bucket
        echo "/usr/bin/time -v aws s3 --profile $awsprofile sync --request-payer requester $s3path/ $pTB/$suffix/ > $logd/s3copy 2>&1" >> copy_finished_pump_runs2opendata.jobs
        echo "$f|$s3path|$suffix" >> $LOG_DIR/all_copies.${date}
    done
    /usr/bin/time -v parallel -j$threads < copy_finished_pump_runs2opendata.jobs > copy_finished_pump_runs2opendata.jobs.run${threads} 2>&1
    echo -n "" > $LOG_DIR/bad_copies.${date}
    echo -n "" > $LOG_DIR/all_removals.${date}
    echo -n "" > remove_finished_pump_runs2opendata.jobs
    for copyline in `cat $LOG_DIR/all_copies.${date}`; do
        f=$(echo "$copyline" | cut -d'|' -f 1)
        s3path=$(echo "$copyline" | cut -d'|' -f 2)
        suffix=$(echo "$copyline" | cut -d'|' -f 3)
        logd=$LOG_DIR/$suffix
        set +eo pipefail
        good=$(fgrep "Exit status: 0" $logd/s3copy)
        set -eo pipefail
        if [[ -n $good ]]; then
            /usr/bin/time -v aws s3 --profile $awsprofile cp --request-payer requester $pDONES/$f $pTB/DONES/ > $logd/s3copy.done 2>&1
            #then remove path from source bucket
            echo "/usr/bin/time -v aws s3 --profile $awsprofile rm --request-payer requester --recursive $s3path/ > $logd/s3rm 2>&1" >> remove_finished_pump_runs2opendata.jobs
        echo "$f|$s3path|$suffix" >> $LOG_DIR/all_removals.${date}
        else
            echo "$copyline" >> $LOG_DIR/bad_copies.${date}
        fi
    done
    /usr/bin/time -v parallel -j$threads < remove_finished_pump_runs2opendata.jobs > remove_finished_pump_runs2opendata.jobs${threads} 2>&1
    echo -n "" > $LOG_DIR/bad_removals.${date}
    for copyline in `cat $LOG_DIR/all_removals.${date}`; do
        f=$(echo "$copyline" | cut -d'|' -f 1)
        s3path=$(echo "$copyline" | cut -d'|' -f 2)
        suffix=$(echo "$copyline" | cut -d'|' -f 3)
        logd=$LOG_DIR/$suffix
        set +eo pipefail
        good=$(fgrep "Exit status: 0" $logd/s3rm)
        set -eo pipefail
        if [[ -n $good ]]; then
            /usr/bin/time -v aws s3 --profile $awsprofile rm --request-payer requester $pDONES/$f > $logd/s3rm.done 2>&1
        else
            echo "$copyline" >> $LOG_DIR/bad_removals.${date}
        fi
    done
    sleep $SLEEP_TIME
done
