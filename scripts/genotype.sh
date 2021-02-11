#!/usr/bin/env bash
#basic genotyping from RNA-seq BAM for sample QC purposes
set -exo pipefail

BAM=$1
SNPS=$2
GENOME_FA=$3

samtools mpileup -l $SNPS -AB -q0 -Q13 -d1000000 -uf $GENOME_FA ${BAM} -o ${BAM}.snps_temp
#~/time $dir/bcftools mpileup --threads 4 -R $snps_file -AB -q0 -Q13 -d1000000 -f $genome_fa -O v -o ${BAM}.snps_temp ${BAM}
bcftools call -mv -Oz ${BAM}.snps_temp > ${BAM}.genotyped.tsv.gz
rm -f ${BAM}.snps_temp
