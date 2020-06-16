#!/usr/bin/env bash

# Author: Chris Wilks & Ben Langmead
#  Email: broadsword@gmail.com
#   Date: 6/16/2020

# Designed to run the aggregator step (recount-unify)
# in the Monorail workflow.  This aggregates across
# multiple samples in a study, project, or compilation

set -ex

#assumes we're in a container where the project-specific
#working directory has been bind mounted and we've entered it

# Ensure directories/paths

#path to symbolic links of all input files/directory structure
#to recount-unify, e.g. ./links
#this needs to be created/populated ahead of running this script
echo "Input dir: ${INPUT_DIR}"
test -n "${INPUT_DIR}"
test -d "${INPUT_DIR}"

#e.g. /path/to/annotated_junctions.tsv.gz
echo "Annotated JXs Path: ${ANNOTATED_JXS}"
test -n "${ANNOTATED_JXS}"
test -s "${ANNOTATED_JXS}"

#e.g. /path/to/exons.bed.w_header.gz
echo "Disjoint Exons BED Path: ${EXON_COORDINATES_BED}"
test -n "${EXON_COORDINATES_BED}"
test -s "${EXON_COORDINATES_BED}"

#e.g. /path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed
echo "Gene rejoin mapping Path: ${GENE_REJOIN_MAPPING}"
test -n "${GENE_REJOIN_MAPPING}"
test -s "${GENE_REJOIN_MAPPING}"

#e.g. /path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed
echo "Gene annotation mapping Path: ${GENE_ANNOTATION_MAPPING}"
test -n "${GENE_ANNOTATION_MAPPING}"
test -s "${GENE_ANNOTATION_MAPPING}"

#e.g. /path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed
echo "Exon rejoin mapping Path: ${EXON_REJOIN_MAPPING}"
test -n "${EXON_REJOIN_MAPPING}"
test -s "${EXON_REJOIN_MAPPING}"

#e.g. /path/to/hg38.recount_pump.fa.new_sizes
echo "Reference FASTA sizes Path: ${REF_SIZES}"
test -n "${REF_SIZES}"
test -s "${REF_SIZES}"

#e.g. /path/to/hg38.recount_pump.fa
echo "Reference FASTA Path: ${REF_FASTA}"
test -n "${REF_FASTA}"
test -s "${REF_FASTA}"

#e.g. G026,G029,R109,ERCC,SIRV,F006
echo "Annotation List: ${LIST_OF_ANNOTATIONS}"
test -n "${LIST_OF_ANNOTATIONS}"

#e.g. 40
echo "CPUs: ${RECOUNT_CPUS}"
test -n "${RECOUNT_CPUS}"

num_samples=$(cat ${SAMPLE_ID_MANIFEST} | wc -l) 

#a config.json file could be provided in running directory
CONFIGFILE=""
if [[ -f "config.json" ]] ; then
    CONFIGFILE="--configfile config.json"
fi
snakemake \
    --snakefile /Snakefile \
    ${CONFIGFILE}
    -j "${RECOUNT_CPUS}" \
    --stats recount-unify.stats.json
    -p
    --config \ 
        input="${INPUT_DIR}" \
        staging=staging \
        compilation_id="${PROJECT_ID}" \
        sample_ids_file="${SAMPLE_ID_MANIFEST}" \
        num_samples=${num_samples} \
        existing_sums="${EXON_COORDINATES_BED}" \
        gene_rejoin_mapping="${GENE_REJOIN_MAPPING}" \
        gene_mapping_final="${GENE_ANNOTATION_MAPPING}" \
        exon_rejoin_mapping="${EXON_REJOIN_MAPPING}" \
        annotated_sjs="${ANNOTATED_JXS}" \
        ref_sizes="${REF_SIZES}" \
        ref_fasta="${REF_FASTA}" \
        annotation_list="${LIST_OF_ANNOTATIONS}" \
        2>&1 | tee recount-unify.output.txt
echo SUCCESS
