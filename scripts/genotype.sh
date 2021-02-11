#!/usr/bin/env bash
#basic genotyping from RNA-seq BAM for sample QC purposes
set -exo pipefail

BAM=$1
SNPS=$2
GENOME_FA=$3
THREADS=$4

if [[ ! -e ${BAM}.bai ]]; then
    samtools index -@ $THREADS $BAM
fi

#samtools mpileup -l $SNPS -AB -q0 -Q13 -d1000000 -uf $GENOME_FA ${BAM} -o ${BAM}.snps_temp
bcftools mpileup --threads $THREADS -R $SNPS -AB -q0 -Q13 -d1000000 -f $GENOME_FA -O v -o ${BAM}.snps_temp ${BAM}
bcftools call -mv -Oz ${BAM}.snps_temp > ${BAM}.genotyped.tsv.gz
rm -f ${BAM}.snps_temp
