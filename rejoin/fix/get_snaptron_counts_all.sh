#!/usr/bin/env bash
set -eo pipefail
dir=$(dirname $0)
threads=30

study_listF=$1

sdir="DONES/"$(echo "$study_listF" | sed 's#studies##')

pushd $sdir

echo -n "" > snaptron_full.jobs
for s in `cat ../../$study_listF`; do
    echo "/bin/bash $dir/convert2snaptron_counts_per_study.sh $s > ${s}.snaptron_full.run 2>&1" >> snaptron_full.jobs
done

echo "pushd $sdir ; /usr/bin/time -v parallel -j${threads} < snaptron_full.jobs > snaptron_full.jobs.run 2>&1 ; wc -l */G026.snaptron | fgrep "63856" | wc -l ; popd"
#wc -l */G026.snaptron | fgrep "3326" | wc -l
#wc -l */G026.snaptron | fgrep -v "3326"
