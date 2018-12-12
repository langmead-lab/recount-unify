#!/usr/bin/env bash
#script to download all annotation sources for junctions
#and rips them to a non-redundant file
#also requires pypy
set -o pipefail -o nounset -o errexit 

date=`date +%b_%d_%Y`
echo $date
annotations_file="${date}_m38_annotations.tar.gz"
echo $annotations_file

mkdir mouse
cd mouse

mkdir -p anno/m38
cd anno/m38
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/knownGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M9/gencode.vM9.basic.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.basic.annotation.gtf.gz
cd ../../
tar -cvzf ${annotations_file} anno

pypy ../rip_annotated_junctions.py --extract-script-dir ../ --annotations ${annotations_file} --org m38 --org-sizes ../m38.sizes 2> m38.annotated_junctions_stats.txt | sort -k1,1 -k2,2n -k3,3n | gzip > annotated_junctions.im38.${date}.tsv.gz
