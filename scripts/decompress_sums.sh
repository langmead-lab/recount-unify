#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit

manifest=$1
output_sentinel=$2
#run one decompress job at a time, but run it for all of the files listed in the manifest
cut -f 3 $manifest | xargs -n 1 -P 1 -I{} sh -c 'zstd -cd $1 | cut -f 4 > ${1}.unc' -- {}
cut -f 3 $manifest | perl -ne 'chomp; $f=$_.".unc"; @s=stat($f); if($s[7]==0) { `cp blank_exon_sums $f`;}'
#when done, write a "done" sentinel
touch $output_sentinel
