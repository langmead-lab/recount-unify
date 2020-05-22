#!/usr/bin/env bash
#script to download all annotation sources for junctions
#and rips them to a non-redundant file
#uses both hg19 and hg38 junctions so needs to have liftover
#also requires pypy
set -o pipefail -o nounset -o errexit 

date=`date +%b_%d_%Y`
echo $date
annotations_file="${date}_hg38_annotations.tar.gz"
echo $annotations_file

LIFTOVER=$1
CHAIN=$2

mkdir human
cd human

mkdir -p anno/hg38
cd anno/hg38
#includes the same as intropolis but with some additions (e.g. CHESS, extra gencodes)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/mgcGenes.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/lincRNAsTranscripts.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/sibGene.txt.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
for i in {20..33}; do
    curl "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${i}/gencode.v${i}.annotation.gtf.gz" > gencode.v${i}.annotation.gtf.gz
    sleep 60
done
#Chess comes in GFF which needs to be slightly modified to be a GTF which extract_splice_sites.py can read
curl "http://ccb.jhu.edu/chess/data/chess2.2_assembly.gff.gz" | zcat | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if($f=~/^#/ || $f[2] ne "exon") { print "$f\n"; next; } $f[8]=~/Parent=((CHS\.\d+)(\.\d+)?)/; $transcript_id=$1; $gene_id=$2; $f[8]="gene_id \"$gene_id\"; transcript_id \"$transcript_id\";"; print "".join("\t",@f)."\n";' | gzip > chess2.2_assembly.gtf.gz
cd ../../
mkdir -p anno/hg19
cd anno/hg19
#exact same list as intropolis
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/mgcGenes.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/lincRNAsTranscripts.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/sibGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/acembly.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
curl "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_3c/gencode.v3c.annotation.GRCh37.gtf.gz" > gencode.v3.annotation.gtf.gz
sleep 60
curl "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_4/gencode.v4.annotation.GRCh37.gtf.gz" > gencode.v4.annotation.gtf.gz
for i in {5..19}; do
    sleep 60
    curl "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${i}/gencode.v${i}.annotation.gtf.gz" > gencode.v${i}.annotation.gtf.gz
done
cd ../../

tar -cvzf ${annotations_file} anno

pypy ../rip_annotated_junctions.py --extract-script-dir ../ --annotations ${annotations_file} --liftover $LIFTOVER --chain $CHAIN --unmapped unmapped_hg19.bed --org hg38 --org-sizes ../hg38.sizes 2> hg38.annotated_junctions_stats.txt | sort -k1,1 -k2,2n -k3,3n | gzip > annotated_junctions.hg38.${date}.tsv.gz
