#!/usr/bin/env bash
#smoke tests of counts of number of studies, rows (genes/exons), and columns (samples) across gene and exon counts files
set -exo pipefail

dir=$(dirname $0)

#env vars should be already defined:
#SAMPLE_ID_MANIFEST
#num_studies
#LIST_OF_ANNOTATIONS
#GENE_EXON_ANNOTATION_ROW_COUNTS_FILE

#need to do smoke tests for proper outputs based on:
#get number of samples per-study
export LC_ALL=C
cut -f 1 $SAMPLE_ID_MANIFEST | sort | uniq -c | tr -s " " \\t | cut -f 2,3 > ${SAMPLE_ID_MANIFEST}.num_samples_per_study.tsv

num_files=$(echo "$LIST_OF_ANNOTATIONS" | tr "," \\n | wc -l)
#2) number of files * number of studies
num_expected=$(( num_studies * num_files ))
num_gene_files=$(find gene_sums_per_study -name "*.gz" -size +0c | wc -l)
if [[ $num_expected -ne $num_gene_files ]]; then
    echo "FAILURE running unify, unexpected # of gene sum files: $num_gene_files vs. $num_expected (expected)"
    exit -1
fi
num_exon_files=$(find exon_sums_per_study -name "*.gz" -size +0c | wc -l)
if [[ $num_expected -ne $num_exon_files ]]; then
    echo "FAILURE running unify, unexpected # of exon sum files: $num_exon_files vs. $num_expected (expected)"
    exit -1
fi

#check 1) row and 2) column (sample) counts per study file for genes
for f in `find gene_sums_per_study -name "*.gz" -size +0c` ; do pcat $f | fgrep -v "##" | perl $dir/check_unifier_outputs.pl $f ${SAMPLE_ID_MANIFEST}.num_samples_per_study.tsv gene ${GENE_EXON_ANNOTATION_ROW_COUNTS_FILE} ; done > gene_checks.tsv

num_passed=$(fgrep "NUM_ROWS_COLUMNS_AND_NO_BLANKS_CHECKS_PASSED" gene_checks.tsv | wc -l)
if [[ $num_passed -ne $num_gene_files ]]; then
    cat <(echo "FAILURE running unify, failed gene sums row/column checks:") gene_checks.tsv
    exit -1
fi

#check 1) row and 2) column (sample) counts per study file for exons
for f in `find exon_sums_per_study -name "*.gz" -size +0c` ; do pcat $f | fgrep -v "##" | perl $dir/check_unifier_outputs.pl $f ${SAMPLE_ID_MANIFEST}.num_samples_per_study.tsv exon ${GENE_EXON_ANNOTATION_ROW_COUNTS_FILE} ; done > exon_checks.tsv

num_passed=$(fgrep "NUM_ROWS_COLUMNS_AND_NO_BLANKS_CHECKS_PASSED" exon_checks.tsv | wc -l)
if [[ $num_passed -ne $num_exon_files ]]; then
    cat <(echo "FAILURE running unify, failed exon sums row/column checks:") exon_checks.tsv
    exit -1
fi
