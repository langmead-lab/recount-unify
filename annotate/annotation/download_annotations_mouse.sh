#!/usr/bin/env bash
#script to download all annotation sources for junctions
#and rips them to a non-redundant file
#also requires pypy
set -o pipefail -o nounset -o errexit 

date=`date +%b_%d_%Y`
echo $date
annotations_file="${date}_m38_annotations.tar.gz"
echo $annotations_file

LIFTOVER=$1
CHAIN=$2

mkdir mouse
cd mouse

mkdir -p anno/m38
cd anno/m38
wget --retry-connrefused ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72311/suppl/GSE72311%5Fmerged%5Ftranscript%5Fassembly%2Egtf%2Egz -O GSE72311_lncrna_pre.gtf.gz
zcat GSE72311_lncrna_pre.gtf.gz | perl -ne '$f=$_; if($f=~/^((\d+)|(MT)|X|Y)\t/) { $f=~s/^MT\t/M\t/; print "chr$f"; next; } print "$f";' | gzip > GSE72311_lncrna.gtf.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ccdsGene.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/knownGene.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/mgcGenes.txt.gz

for i in {2..24}; do
    #wget --retry-connrefused ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${i}/gencode.vM${i}.annotation.gtf.gz
    curl "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${i}/gencode.vM${i}.annotation.gtf.gz" > gencode.vM${i}.annotation.gtf.gz
    sleep 60
done
cd ../../

mkdir -p anno/m37
cd anno/m37
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/ccdsGene.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/knownGene.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/mgcGenes.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/acembly.txt.gz
wget --retry-connrefused http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/vegaGene.txt.gz
cd ../../

tar -cvzf ${annotations_file} anno

pypy ../rip_annotated_junctions.py --extract-script-dir ../ --annotations ${annotations_file} --liftover $LIFTOVER --chain $CHAIN --unmapped unmapped_m37.bed --org m38 --org-sizes ../m38.sizes 2> m38.annotated_junctions_stats.txt | sort -k1,1 -k2,2n -k3,3n | gzip > annotated_junctions.m38.${date}.tsv.gz
