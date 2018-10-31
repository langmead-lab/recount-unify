#!/usr/bin/env bash
#script to download all annotation sources for junctions
#and rips them to a non-redundant file
#uses both hg19 and hg38 junctions so needs to have liftover
#also requires pypy
set -o pipefail -o nounset -o errexit 

date=`date +%b_%d_%Y`
echo $date
annotations_file="${date}_annotations.tar.gz"
echo $annotations_file

LIFTOVER=$1
CHAIN=$2

mkdir -p anno/hg38
cd anno/hg38
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/mgcGenes.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/lincRNAsTranscripts.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/sibGene.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
#Chess comes in GFF which needs to be slightly modified to be a GTF which extract_splice_sites.py can read
curl "http://ccb.jhu.edu/chess/data/chess2.1_assembly.gff.gz" | zcat | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if($f=~/^#/ || $f[2] ne "exon") { print "$f\n"; next; } $f[8]=~/Parent=((CHS\.\d+)(\.\d+)?)/; $transcript_id=$1; $gene_id=$2; $f[8]="gene_id \"$gene_id\"; transcript_id \"$transcript_id\";"; print "".join("\t",@f)."\n";' | gzip > chess2.1_assembly.gtf.gz
cd ../../
mkdir -p anno/hg19
cd anno/hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/mgcGenes.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/lincRNAsTranscripts.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/sibGene.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/acembly.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz
cd ../../

tar -cvzf ${annotations_file} anno

pypy hg38_rip_annotated_junctions.py --hisat2-dir ./ --annotations ${annotations_file} --liftover $LIFTOVER --chain $CHAIN --unmapped unmapped_hg19.bed 2> annotated_junctions_stats.txt | sort -k1,1 -k2,2n -k3,3n | gzip > annotated_junctions.tsv.gz
