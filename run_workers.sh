#!/usr/bin/env bash
#Main multi-Unifier runner script:
#1) formats and mounts SSDs
#2) downloads Monorail Unifier references to SSD
#3) stars N Unifier workers in parallel
#Assumes the following ENV vars are set:
#a) REF (hg38 or grcm38, required)
#b) NUM_WORKERS (number of concurrent workers to start, default: 16)
#c) NUM_CORES (maximum number of CPUs to use per worker, default: 8)
#d) SSD_MIN_SIZE (minimum size of local SSDs, default: 600GBs)
#OPTIONAL:
#KEEP_RUNNING=1 means this script will keep running and polling the queue indefinitely
set -exo pipefail
dir=$(dirname $0)
if [[ -n $DEBUG ]]; then
    sleep 1d
fi
DEFAULT_NUM_WORKERS=16
DEFAULT_SSD_MIN_SIZE=600000000000
DEFAULT_NUM_CORES=8
if [[ -z $REF ]]; then
    echo "no REF set, terminating early!"
    exit -1
fi

#check for # of processors to determine instance type being run on:
num_procs=$(fgrep processor /proc/cpuinfo | wc -l)
if [[ $num_procs -eq 32 ]]; then
    DEFAULT_NUM_WORKERS=8
fi
if [[ $num_procs -eq 48 ]]; then
    DEFAULT_NUM_WORKERS=12
fi

if [[ -z $NUM_WORKERS ]]; then
    export NUM_WORKERS=$DEFAULT_NUM_WORKERS
fi
if [[ -z $NUM_CORES ]]; then
    export NUM_CORES=$DEFAULT_NUM_CORES
fi
if [[ -z $SSD_MIN_SIZE ]]; then
    export SSD_MIN_SIZE=$DEFAULT_SSD_MIN_SIZE
fi

user=$(whoami)
set +eo pipefail
df=$(df | fgrep "/md1" | fgrep "/work1" | wc -l)
set -eo pipefail
#check for local SSDs, creates file local_disks.txt
if [[ $df -eq 0 ]]; then
    #1) format and mount SSDs (but skip root)
    set +eo pipefail
    MAKE_1_FS=1
    sudo /usr/bin/time -v /bin/bash -x $dir/check_and_create_fs_for_ephemeral_disks.sh $MAKE_1_FS
    set -eo pipefail

    num_ssds=$(cat local_disks.txt | wc -l)
    if [[ $num_ssds -eq 0 ]]; then
        export NO_SSD=1
        sudo mkdir /work1
        sudo chown ubuntu /work1
        sudo chmod u+rwx /work1
    fi
fi

#2) download Unifier references to SSD
if [[ ! -d /work1/ref/${REF}_unify ]]; then
    mkdir -p /work1/ref
    pushd /work1/ref
    #/usr/bin/time -v $dir/../get_unify_refs.sh $REF > get_unify_refs.sh.${REF}.run 2>&1
    /usr/bin/time -v $dir/get_unify_refs.additional.sh $REF > get_unify_refs.additional.sh.${REF}.run 2>&1
    popd
fi 

#3) start N concurrent Unifier workers
export REF_DIR=/work1/ref
pushd /work1
mkdir -p /work1/runs
if [[ $num_ssds -gt 1 ]]; then
    mkdir -p /work2/runs
fi
echo -n "" > worker.jobs
idx=2
for i in $( seq 1 $NUM_WORKERS ); do 
    #alternate SSD to write job on
    if [[ $idx -eq 1 && $num_ssds -gt 1 ]]; then
        idx=2
    else
        idx=1
    fi
    echo "/usr/bin/time -v /bin/bash -x $dir/worker.sh /work${idx} > /work1/runs/w${i}.run 2>&1" >> worker.jobs
done

#setup main SRA metadata file that's already pre-compiled and on S3
#TODO update to later date in 2024
export SRA_METADATA_PRECOMPILED_S3_URL="s3://monorail-batch/all_runs.tsv.2007_to_20230422.nodups.nomissing_studies.tsv.gz"
SRA_METADATA_BASE=$(basename $SRA_METADATA_PRECOMPILED_S3_URL)
export SRA_METADATA=$(echo "/work1/$SRA_METADATA_BASE" | sed 's#\.gz$##')
aws s3 cp $SRA_METADATA_PRECOMPILED_S3_URL - | pcat > $SRA_METADATA
export SKIP_FETCHING_SRA_METADATA=1
cat worker.jobs
/usr/bin/time -v parallel -j${NUM_WORKERS} < worker.jobs > worker.jobs.run 2>&1
