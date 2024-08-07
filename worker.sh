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
#g) S3_OUTPUT (where to upload the unifier results to on S3, default: s3://monorail-batch/unifier_outputs/)
#OPTIONAL:
#KEEP_RUNNING=1 means this script will keep running and polling the queue indefinitely
set -exo pipefail
dir=$(dirname $0)

#filesystem to write temporary output to
#e.g. /work2 (defual: /work1)
fs=$1

if [[ -z $fs ]]; then
    export fs="/work1"
fi

if [[ -z $REF ]]; then
    echo "no REF set, terminating early!"
    exit -1
fi
if [[ -z $Q ]]; then
    export Q="https://sqs.us-east-1.amazonaws.com/315553526860/monorail_batch_unify"
fi
export REGION=$(echo "$Q" | cut -d'.' -f 2)
if [[ -z $DOCKER_IMAGE ]]; then
    export DOCKER_IMAGE="315553526860.dkr.ecr.us-east-2.amazonaws.com/recount-unify-aarch64:1.1.5"
fi
if [[ -z $REF_DIR ]]; then
    export REF_DIR=/work1/ref
fi
if [[ -z $NUM_CORES ]]; then
    export NUM_CORES=8
fi
if [[ -z $OUTPUT_DIR_GLOBAL ]]; then
    export OUTPUT_DIR_GLOBAL="$fs/unifier"
fi
if [[ -z $S3_OUTPUT ]]; then
    #test destination bucket
    #export S3_OUTPUT="s3://monorail-batch/unifier-output"
    #production run to AWS Open Data live bucket!:
    #export S3_OUTPUT="s3://recount-opendata/recount3/release/human/data_sources/sra"
    #export S3_OUTPUT="s3://recount-opendata/recount3/release/human/data_sources/sra"
    org0="human"
    if [[ "$REF" == "grcm38" ]]; then
        org0="mouse"
    fi
    #export S3_OUTPUT="s3://recount-opendata/recount3expansion/unifier/$org0"
    export S3_OUTPUT="s3://monorail-batch/unifier-output/$org0"
fi
export S3_UNIFIER_DONES="s3://monorail-batch/UNIFIER_DONES"

#1) check for new studies on the queue
msg_json=$(aws sqs receive-message --region $REGION --queue-url $Q)
while [[ -n $msg_json || -n $KEEP_RUNNING ]]; do
    if [[ -n $msg_json ]]; then
        set +eo pipefail
        handle=$(echo "$msg_json" | fgrep '"ReceiptHandle":' | cut -d'"' -f 4)
        #e.g. s3://monorail-batch/pump_outputs/<lo>/<study>
        study_s3=$(echo "$msg_json" | fgrep '"Body":' | cut -d'"' -f 4)
        set -eo pipefail
        if [[ -z $handle || -z $study_s3 ]]; then
            echo "ERROR: didn't find either a handle or a study in SQS message: $msg_json.  skipping"
            aws sqs delete-message --region $REGION --queue-url $queue --receipt-handle $handle
            msg_json=$(aws sqs receive-message --region $REGION --queue-url $Q)
            continue
        fi
        date=$(date +%Y%m%d_%s)
        study=$(basename $study_s3)
        lo=${study: -2}
        export OUTPUT_DIR=$OUTPUT_DIR_GLOBAL/$study.${date}
        rm -rf $OUTPUT_DIR
        mkdir -p $OUTPUT_DIR
        pushd $OUTPUT_DIR
        #TODO: start a unifier job on the study
        #2) download from S3 pump outputs for study
        #/usr/bin/time -v aws s3 cp --recursive $study_s3 $study/
        #2a) get SRA metadata for study from pre-compiled file on S3 to avoid having to query SRA for it per-study
        #fgrep $'\t'"$study"$'\t' $SRA_METADATA 
        #2b)
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION START
        echo -n "" > ${study}.s3dnload.jobs
        echo -n "" > samples4study.tsv.temp
        bucket=$(echo "$study_s3" | cut -d'/' -f 3)
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION CREATING_PUMP_OUTPUT_DOWNLOAD_JOBS
        num_bws2dnload=0 
        for f in `aws s3 ls --no-sign-request --recursive $study_s3 | fgrep "manifest" | tr -s " " $'\t' | cut -f4`; do
            num_bws2dnload=$((num_bws2dnload+1))
            f0=$(basename $f | cut -d'!' -f1)
            study_=$(basename $f | cut -d'!' -f 2)
            echo "	$study_	$f0	" >> samples4study.tsv.temp
            s3path=$(dirname $f)
            sample=$(basename $s3path) 
            #don't need the largest files---bigwigs nor nonrefs---for Unifier
            #but download the # of bigwigs that equal the number of cores to check
            if [[ $num_bws2dnload -le $NUM_CORES ]]; then
                echo "/usr/bin/time -v aws s3 cp --no-sign-request --recursive --exclude \"*.unique.bw\" --exclude \"*.bamcount_nonref.csv.zst\" s3://$bucket/$s3path/ ./$sample/ > ../runs/${sample}.s3dnload 2>&1" >> ${study}.s3dnload.jobs
            else
                echo "/usr/bin/time -v aws s3 cp --no-sign-request --recursive --exclude \"*.unique.bw\" --exclude \"*.all.bw\" --exclude \"*.bamcount_nonref.csv.zst\" s3://$bucket/$s3path/ ./$sample/ > ../runs/${sample}.s3dnload 2>&1" >> ${study}.s3dnload.jobs
            fi
        done
        #echo $'study_id\tsample_id' > samples4study.tsv
        head -1 $SRA_METADATA | cut -f2- > samples4study.tsv
        fgrep -f samples4study.tsv.temp $SRA_METADATA | cut -f 2- | LC_ALL=C sort -u >> samples4study.tsv
        wc_expected=$(fgrep $'\t'"$study_"$'\t' $SRA_METADATA | wc -l)
        wc0=$(cat samples4study.tsv | wc -l)
        wc1=$(cat samples4study.tsv.temp | wc -l)
        #add +1 for the header
        wc1=$((wc1 + 1))
        wc_expected=$((wc_expected + 1))
        if [[ $wc_expected -ne $wc0 || $wc0 -ne $wc1 ]]; then
            echo "number of samples/runs does not match between precompiled SRA metadata and pump run: $wc_expected vs. $wc0 vs. $wc1 from $study_s3, skipping"
            /bin/bash $dir/monorail_unifier_log.sh $study $REGION UNIFIER_METADATA_PULL_FAILED
            msg_json=$(aws sqs receive-message --region $REGION --queue-url $Q)
            continue
        fi
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION METADATA_PULL_GOOD
        rm -f samples4study.tsv.temp
        mkdir -p runs
        if [[ ! -d pump ]]; then
            mkdir -p pump
            pushd pump
            /usr/bin/time -v parallel -j${NUM_CORES} < ../${study}.s3dnload.jobs > ../${study}.s3dnload.jobs.run 2>&1
            popd
        fi
        num_samples=$(cat ${study}.s3dnload.jobs | wc -l)
        num_downloaded=$(fgrep "Exit " runs/*.s3dnload | fgrep 'Exit status: 0' | wc -l)
        if [[ $num_downloaded -ne $num_samples ]]; then
            echo "not all were able to download, skipping study $study"
            /bin/bash $dir/monorail_unifier_log.sh $study $REGION UNIFIER_PUMP_OUTPUT_DNLOAD_FAILED
            msg_json=$(aws sqs receive-message --region $REGION --queue-url $Q)
            continue
        fi
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION PUMP_OUTPUT_DNLOAD_GOOD
        #3) Run Unifier
        #TODO: double check params
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION UNIFY_PROPER_START
        /usr/bin/time -v /bin/bash -x $dir/run_recount_unify_within_container.sh $REF $REF_DIR `pwd`/unifier `pwd`/pump `pwd`/samples4study.tsv $NUM_CORES > run_recount_unify_within_container.sh.run 2>&1
        success=$(egrep -e '^SUCCESS$' run_recount_unify_within_container.sh.run)
        if [[ -z $success ]]; then
            echo "unifier failed for $study, skipping"
            /bin/bash $dir/monorail_unifier_log.sh $study $REGION UNIFIER_PROPER_FAILED
            msg_json=$(aws sqs receive-message --region $REGION --queue-url $Q)
            continue
        fi
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION UNIFIER_PROPER_DONE
        #4) copy unifier results back to S3
        pushd `pwd`/unifier/run_files
        mv all.logs.tar.gz ../
        rm -rf staging_jxs staging input_from_pump links *.pre_existing *.gz blank_exon_sums *.annotation.tsv *.gene_sums.tsv intron_counts_summed.tsv *.RR *.mm *.coords ../junction_counts_per_study_run_files
        popd
        #/usr/bin/time -v aws s3 cp --recursive `pwd`/unifier/ $S3_OUTPUT/$lo/$study.${date}/ > s3upload.run 2>&1
        #UPDATE: push back to target AWS Open Data recount3 bucket release structure (this is production)
        echo -n "" > s3upload.run
        for d0 in gene_sums exon_sums junctions metadata; do
            ds0="$d0"
            if [[ $d == "gene_sums" || $d == "exon_sums" ]]; then
                ds0="${d}_per_study"
            elif [[ $d == "junctions" ]]; then
                ds0="${d}_counts_per_study"
            fi
            /usr/bin/time -v aws s3 cp --recursive `pwd`/unifier/ $S3_OUTPUT/$d0/$lo/$study/ >> s3upload.run 2>&1
        done
        /usr/bin/time -v aws s3 cp `pwd`/unifier/all.logs.tar.gz $S3_OUTPUT/unifier_logs/$lo/$study/ >> s3upload.run 2>&1
        echo "$S3_OUTPUT/$d0/$lo/$study" > ${study}.DONE
        aws s3 cp ${study}.DONE $S3_UNIFIER_DONES/
        /bin/bash $dir/monorail_unifier_log.sh $study $REGION END
        #get next message repeat
        aws sqs delete-message --region $REGION --queue-url $Q --receipt-handle $handle
        popd
        #final cleanup
        rm -rf $OUTPUT_DIR
        #aws s3 rm --recursive $study_s3/
    else
        sleep 60
    fi
    msg_json=$(aws sqs receive-message --region $REGION --queue-url $Q)
done
