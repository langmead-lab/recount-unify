#!/bin/bash

ref=$1

if [[ ! -v $ref ]]; then ref=hg38 ; fi

echo "Downloading ${ref} Junctions..."
aws --profile jhu-langmead s3 cp s3://recount-ref/${ref}_junctions/annotated_junctions.tar.gz ./
tar -zxvf annotated_junctions.tar.gz
