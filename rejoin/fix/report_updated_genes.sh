#!/usr/bin/env bash
#set -eo pipefail

study=$1
src=$2

for annot in G026 G029 R109 F006; do
    num=$(fgrep "<" $study/${annot}.diff | cut -d' ' -f2- | cut -f1 | sort -u | wc -l)
    if [[ $num -gt 0 ]]; then
        #pcat $study/${src}.gene_sums.${study}.${annot}.gz
        echo -n "$study	$annot	$num	,"
        fgrep "<" $study/${annot}.diff | cut -d' ' -f2- | cut -f1 | sort -u | tr $'\n' ","
        echo ""
    fi
done
