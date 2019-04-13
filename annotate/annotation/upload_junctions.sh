#!/bin/bash

ref=$1

if [[ ! -v $ref ]]; then ref=hg38 ; fi

junctions=`ls -t annotated_junctions.${ref}.*.tsv.gz | head -1`

ln -fs $junctions annotated_junctions.tsv.gz
tar -zcvf annotated_junctions.tar.gz annotated_junctions.tsv.gz $junctions

echo "Uploading ${ref} Junctions..."
aws --profile jhu-langmead s3 cp annotated_junctions.tar.gz s3://recount-ref/${ref}_junctions/
