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
if [[ -z $REF ]]; then
    echo "no REF set, terminating early!"
    exit -1
fi
if [[ -z $NUM_WORKERS ]]; then
    export NUM_WORKERS=16
fi
if [[ -z $NUM_CORES ]]; then
    export NUM_CORES=8
fi
if [[ -z $SSD_MIN_SIZE ]]; then
    export SSD_MIN_SIZE=600000000000
fi

df=$(df | fgrep "/work1")
#1) format and mount SSDs (but skip root)
if [[ -z $df ]]; then
    i=1
    for d in `lsblk -b | egrep -e '^nvme' | fgrep -v nvme0n1 | perl -ne '@f=split(/\s+/,$_,-1); next if($f[3] < '$SSD_MIN_SIZE'); print $f[0]."\n";'`; do sudo mkfs -q -t ext4 /dev/$d ; sudo mkdir -p /work${i} ; sudo mount /dev/$d /work${i}/ ; sudo chown -R ubuntu /work${i} ; sudo chmod -R a+rw /work${i} ; i=$((i + 1)); done
fi

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
echo -n "" > worker.jobs
for i in $( seq 1 $NUM_WORKERS ); do echo "/usr/bin/time -v /bin/bash -x $dir/worker.sh > /work1/runs/w${i}.run 2>&1" >> worker.jobs ; done
cat worker.jobs
/usr/bin/time -v parallel -j${NUM_WORKERS} < worker.jobs > worker.jobs.run 2>&1
