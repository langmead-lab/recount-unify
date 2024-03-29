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

set +eo pipefail
df=$(df | fgrep "/work" | wc -l)
set -eo pipefail
#1) format and mount SSDs (but skip root)
user=$(whoami)
i=1
if [[ $df -eq 0 ]]; then
    for d in `lsblk -b | egrep -e '^nvme' | fgrep -v nvme0n1 | perl -ne '@f=split(/\s+/,$_,-1); next if($f[3] < '$SSD_MIN_SIZE'); print $f[0]."\n";'`; do 
        sudo mkfs -q -t ext4 /dev/$d
        sudo mkdir -p /work${i}
        sudo mount /dev/$d /work${i}/
        sudo chown -R $user /work${i}
        sudo chmod -R a+rw /work${i}
        i=$((i + 1))
    done
else
    i=$((df + 1))
fi
num_ssds=$((i - 1))

#2) download Unifier references to SSD
if [[ ! -d /work1/ref/${REF}_unify ]]; then
    mkdir -p /work1/ref
    pushd /work1/ref
    /usr/bin/time -v $dir/../get_unify_refs.sh $REF > get_unify_refs.sh.${REF}.run 2>&1
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
export SRA_METADATA_PRECOMPILED_S3_URL="s3://monorail-batch/all_runs.tsv.2007_to_20230422.nodups.nomissing_studies.tsv.gz"
SRA_METADATA_BASE=$(basename $SRA_METADATA_PRECOMPILED_S3_URL)
export SRA_METADATA=$(echo "/work1/$SRA_METADATA_BASE" | sed 's#\.gz$##')
aws s3 cp $SRA_METADATA_PRECOMPILED_S3_URL - | pcat > $SRA_METADATA
export SKIP_FETCHING_SRA_METADATA=1
cat worker.jobs
/usr/bin/time -v parallel -j${NUM_WORKERS} < worker.jobs > worker.jobs.run 2>&1
