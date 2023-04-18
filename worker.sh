#!/usr/bin/env bash
#Main single-Unifier runner script:
#in a loop:
#1) checks for job on unifier SQS (one study)
#2) dequeues job and gets date_time stamp at time of dequeue
#3) copies down from S3 all results and setups up local FS links
#4) calls workflow.bash.recount_jxns_only (recount3 only) on study
#5) on success, copies results into final S3 location (Open Data) and either removes or moves the pump results (to avoid re-queuing)
#6) regardless of status of 4/5, removes results and temporary files from local SSDs to keep usage down
#7) repeat
#Assumes the following ENV vars are set:
#a) REF (hg38 or grcm38, required)
#b) Q (queue, usually SQS URL) to retrieve jobs from (1 per study to unify, default: monorail_batch on SQS)
#c) DOCKER_IMAGE (docker image of Unifier to use, default: 1.1.1)
#d) REF_DIR (references directory, default: /work1/ref)
#e) NUM_CORES (maximum number of CPUs to use per worker, default: 8)
#f) OUTPUT_DIR_GLOBAL (where to write the unifier outputs temporarily, before copying back to S3, default: /work1/unifier/$study)
#g) S3_OUTPUT (where to upload the unifier results to on S3, default: s3://neuro-datalake-ds/research/genomics/raw_data/monorail_temp_input/unifier_outputs/)
set -exo pipefail
dir=$(dirname $0)

if [[ -z $REF ]]; then
    echo "no REF set, terminating early!"
    exit -1
fi
if [[ -z $Q ]]; then
    export Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_unify"
fi
if [[ -z $DOCKER_IMAGE ]]; then
    export DOCKER_IMAGE="315553526860.dkr.ecr.us-east-1.amazonaws.com/recount-unify-aarch64:1.1.1"
fi
if [[ -z $REF_DIR ]]; then
    export REF_DIR=/work1/ref
fi
if [[ -z $NUM_CORES ]]; then
    export NUM_CORES=8
fi
if [[ -z $OUTPUT_DIR_GLOBAL ]]; then
    export OUTPUT_DIR_GLOBAL=/work2/unifier
fi
if [[ -z $S3_OUTPUT ]]; then
    export S3_OUTPUT="s3://neuro-datalake-ds/research/genomics/raw_data/monorail_temp_input/unifier_outputs"
fi

#1) check for new studies on the queue
msg_json=$(aws sqs receive-message --queue-url $Q)
while [[ -n $msg_json ]]; do
    set +eo pipefail
    handle=$(echo "$msg_json" | fgrep '"ReceiptHandle":' | cut -d'"' -f 4)
    #e.g. s3://neuro-datalake-ds/research/genomics/raw_data/monorail_temp_input/pump_outputs/10/SRP277410
    study_s3=$(echo "$msg_json" | fgrep '"Body":' | cut -d'"' -f 4)
    set -eo pipefail
    if [[ -z $handle || -z $study_s3 ]]; then
        echo "ERROR: didn't find either a handle or a study in SQS message: $msg_json.  skipping"
        aws sqs delete-message --queue-url $queue --receipt-handle $handle
        msg_json=$(aws sqs receive-message --queue-url $Q)
        continue
    fi
    date=$(date +%Y%m%d_%s)
    study=$(basename $study_s3)
    lo=${study: -2}
    if [[ -z $OUTPUT_DIR ]]; then
        export OUTPUT_DIR=$OUTPUT_DIR_GLOBAL/$study.${date}
        mkdir -p $OUTPUT_DIR
    fi
    pushd $OUTPUT_DIR
    #TODO: start a unifier job on the study
    #2) download from S3 pump outputs for study
    #/usr/bin/time -v aws s3 cp --recursive $study_s3 $study/
    echo -n "" > ${study}.s3dnload.jobs
    echo -n "" > samples4study.tsv.temp
    bucket=$(echo "$study_s3" | cut -d'/' -f 3)
    for f in `aws s3 ls --recursive $study_s3 | fgrep "manifest" | tr -s " " $'\t' | cut -f4`; do
        f0=$(basename $f | cut -d'!' -f1)
        study_=$(basename $f | cut -d'!' -f 2)
        echo "$study_	$f0" >> samples4study.tsv.temp
        s3path=$(dirname $f)
        sample=$(basename $s3path) 
        echo "/usr/bin/time -v aws s3 cp --recursive s3://$bucket/$s3path/ ./$sample/ > ../runs/${sample}.s3dnload 2>&1" >> ${study}.s3dnload.jobs
    done
    echo $'study_id\tsample_id' > samples4study.tsv
    LC_ALL=C sort -u samples4study.tsv.temp >> samples4study.tsv
    rm -f samples4study.tsv.temp
    mkdir -p runs
    mkdir -p unify
    if [[ ! -d pump ]]; then
        mkdir pump
        pushd pump
        /usr/bin/time -v parallel -j${NUM_CORES} < ../${study}.s3dnload.jobs > ../${study}.s3dnload.jobs.run 2>&1
        popd
    fi
    num_samples=$(cat ${study}.s3dnload.jobs | wc -l)
    num_downloaded=$(fgrep "Exit " runs/*.s3dnload | fgrep 'Exit status: 0' | wc -l)
    if [[ $num_downloaded -ne $num_samples ]]; then
        echo "not all were able to download, skipping study $study"
        msg_json=$(aws sqs receive-message --queue-url $Q)
        continue
    fi
    #3) Run Unifier
    #TODO: double check params
    /usr/bin/time -v /bin/bash -x $dir/run_recount_unify_within_container.sh $REF $REF_DIR `pwd`/unifier `pwd`/pump `pwd`/samples4study.tsv $NUM_CORES > run_recount_unify_within_container.sh.run 2>&1
    success=$(egrep -e '^SUCCESS$' run_recount_unify_within_container.sh.run)
    if [[ -z $success ]]; then
        echo "unifier failed for $study, skipping"
        msg_json=$(aws sqs receive-message --queue-url $Q)
        continue
    fi
    #4) copy unifier results back to S3
    pushd `pwd`/unifier/run_files
    mv all.logs.tar.gz ../
    rm -rf staging_jxs staging input_from_pump links *.pre_existing *.gz blank_exon_sums *.annotation.tsv *.gene_sums.tsv intron_counts_summed.tsv *.RR *.mm *.coords ../junction_counts_per_study_run_files
    popd
    /usr/bin/time -v aws s3 cp --recursive `pwd`/unifier/ $S3_OUTPUT/$lo/$study.${date}/ > s3upload.run 2>&1
    #get next message repeat
    aws sqs delete-message --queue-url $Q --receipt-handle $handle
    msg_json=$(aws sqs receive-message --queue-url $Q)
done
