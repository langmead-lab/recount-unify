#!/usr/bin/env bash
set -exo pipefail 
annot="G026"

study_listF=$1

sdir="DONES/"$(echo "$study_listF" | sed 's#studies##')

pushd $sdir

pcmd="paste "
for s in `cat ../../$study_listF`; do
    #pcmd="${pcmd} $s/${annot}.snaptron"
    pcmd="${pcmd} $s/${annot}.snaptron_full"
done 
$pcmd | sed 's#\t##g' > ../../${study_listF}.snaptron.pasted.1
#paste /datascope/recount03/rejoin_fix/disjoint2exons2genes.fix.sorted.bed.G026_genes.tabs.sorted ../../${study_listF}.snaptron.pasted.1 > ../../${study_listF}.snaptron.pasted

popd
