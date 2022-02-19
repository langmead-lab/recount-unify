#!/usr/bin/env bash
set -exo pipefail
dir=$(dirname $0)
threads=30

study_listF=$1

sdir="DONES/"$(echo "$study_listF" | sed 's#studies##')

pushd $sdir

echo -n "" > snaptron.jobs
for s in `cat ../../$study_listF`; do
    echo "/bin/bash $dir/get_updated_snaptron_counts_per_study.sh $s > ${s}.snaptron.run 2>&1" >> snaptron.jobs
done

/usr/bin/time -v parallel -j${threads} < snaptron.jobs > snaptron.jobs.run 2>&1
wc -l */G026.snaptron | fgrep "3326" | wc -l
wc -l */G026.snaptron | fgrep -v "3326"
