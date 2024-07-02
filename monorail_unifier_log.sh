#!/usr/bin/env bash
#set -exo pipefail
export MD_IP=169.254.169.254
export LOG_GROUP="monorail-unifier"
export LOG_STREAM="unifier_runs"
study=$1
#checkpoint OR path2fastqfiles,
#checkpoint is one of:
#+) START
#+) CREATING_PUMP_OUTPUT_DOWNLOAD_JOBS
#+) METADATA_PULL_GOOD
#+) (or) UNIFIER_METADATA_PULL_FAILED
#+) PUMP_OUTPUT_DNLOAD_GOOD
#+) (or) UNIFIER_PUMP_OUTPUT_DNLOAD_FAILED
#+) UNIFY_PROPER_START
#+) UNIFIER_PROPER_DONE (before copy back to S3)
#+) (or) UNIFIER_PROPER_FAILED
#+) END (after copy back to S3)
mode=$2

#need to establish a unique ID for this instance of processing this RUN accession through unifier:
#date.run_acc.study_acc.node_ip
DATE=$(date +"%Y%m%dT%H%M%S%z")
if [[ ! -e /dev/shm/INSTANCE_INFO ]]; then
    IP=$(curl http://${MD_IP}/latest/meta-data/local-ipv4)
    itype=$(curl http://${MD_IP}/latest/meta-data/instance-type)
    #IP=$(echo "$IP" | sed 's#\.#-#g')
else
    instance_info=$(cat /dev/shm/INSTANCE_INFO)
    IP=$(echo "$instance_info" | cut -d';' -f1)
    itype=$(echo "$instance_info" | cut -d';' -f2)
fi
JOB_ID="${study}"

#collect additional stats from machine
#disk usage
df=$(df -h | tr -s " " $'|' | cut -d'|' -f 4,5,6 | tail -n+2 | fgrep -v "|/snap/" | fgrep -v "|/run" | fgrep -v "|/sys" | fgrep -v "|/dev" | tr $'\n' ";" | sed 's#;$#\n#')
ncores=$(fgrep -i processor /proc/cpuinfo  | wc -l) 
ram=$(head -3 /proc/meminfo | tr $'\n' ";" | sed 's# ##g' | sed 's#;$#\n#')
load_avg=$(top -b -n  1  | awk '/load average/ { printf "%s\n", $12 }' | sed 's#,##')

#now log:
entry="${DATE};${mode};${JOB_ID};${IP};${itype};${ncores};${load_avg};${ram};${df}"
echo "$entry"
d2=$(($(date +%s)*1000))
aws logs put-log-events --log-group-name $LOG_GROUP --log-stream-name $LOG_STREAM --log-events "timestamp=${d2},message=${entry}"
