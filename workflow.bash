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
echo "Working dir: ${WORKING_DIR}"
test -n "${WORKING_DIR}"
test -d "${WORKING_DIR}"

pushd ${WORKING_DIR}

#needed for blank sums file copy
echo "Ref dir: ${REF_DIR}"
test -n "${REF_DIR}"
test -d "${REF_DIR}"

cp ${REF_DIR}/blank_exon_sums ${WORKING_DIR}/

#original root directory of recount-pump output
#OR
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
echo "Disjoint Exons BED Path (w/ header): ${EXON_COORDINATES_BED}"
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

#e.g. srav3h.exon_counts.bitmasks.tsv
echo "Exon Bitmasks Path: ${EXON_BITMASKS}"
test -n "${EXON_BITMASKS}"
test -s "${EXON_BITMASKS}"

#e.g. srav3h.exon_counts.bitmasks.coords
echo "Exon Bitmask Coordinates Path: ${EXON_BITMASK_COORDS}"
test -n "${EXON_BITMASK_COORDS}"
test -s "${EXON_BITMASK_COORDS}"

#e.g. G026,G029,R109,F006,ERCC,SIRV
echo "Annotation List: ${LIST_OF_ANNOTATIONS}"
test -n "${LIST_OF_ANNOTATIONS}"

if [[ -z $PROJECT_ID ]]; then
    PROJECT_ID=0
fi

if [[ -z $PROJECT_SHORT_NAME ]]; then
    PROJECT_SHORT_NAME="sra"
fi

echo "Sample ID manifest: ${SAMPLE_ID_MANIFEST}"
test -n "${SAMPLE_ID_MANIFEST}"
test -e "${SAMPLE_ID_MANIFEST}"

#check passed in sample ID file for number of fields (no header)
num_cols=$(head -1 ${SAMPLE_ID_MANIFEST} | tr \\t \\n | wc -l)
#not getting the numeric IDs (3rd column)
#so we'll need to generate them
if [[ $num_cols -ne 3 ]]; then
    if [[ $num_cols -ne 2 ]]; then
        echo "${SAMPLE_ID_MANIFEST} needs to have either 2 or 3 columns, exiting"
        exit -1
    fi
    pushd /recount-unify/sample_ids
    /usr/bin/python2.7 assign_compilation_ids.py --accessions-file $SAMPLE_ID_MANIFEST --acc-col 1 --compilation-code $PROJECT_ID --no-header > ${WORKING_DIR}/ids.tsv 2> ${WORKING_DIR}/assign_compilation_ids.py.errs
    export SAMPLE_ID_MANIFEST=${WORKING_DIR}/ids.tsv
    popd
fi

/bin/bash -x /recount-unify/scripts/create_directory_hierarchy_for_one_study.sh $SAMPLE_ID_MANIFEST $INPUT_DIR ${WORKING_DIR}/intermediate_links > setup_intermediate_links.run 2>&1

#need to create the directory structure we expect from the basic output of
#recount-pump on one study (assumed)
/bin/bash -x /recount-unify/scripts/find_done.sh ${WORKING_DIR}/intermediate_links links "*_att" > setup_links.run 2>&1

#e.g. 40
echo "CPUs: ${RECOUNT_CPUS}"
test -n "${RECOUNT_CPUS}"

num_samples=$(cat ${SAMPLE_ID_MANIFEST} | wc -l)
num_exons=$(zcat ${EXON_COORDINATES_BED} | tail -n+2 | wc -l)

#a config.json file could be provided in running directory
CONFIGFILE=""
if [[ -f "config.json" ]] ; then
    CONFIGFILE="--configfile config.json"
fi

#first do exon/gene sums Snakemake
snakemake \
    --snakefile /recount-unify/Snakefile \
    ${CONFIGFILE} \
    -j "${RECOUNT_CPUS}" \
    --stats recount-unify.sums.stats.json \
    -p \
    --config \
        input=links \
        staging=staging \
        compilation="${PROJECT_SHORT_NAME}" \
        sample_ids_file="${SAMPLE_ID_MANIFEST}" \
        num_samples=${num_samples} \
        existing_sums="${EXON_COORDINATES_BED}" \
        gene_rejoin_mapping="${GENE_REJOIN_MAPPING}" \
        gene_mapping_final="${GENE_ANNOTATION_MAPPING}" \
        exon_rejoin_mapping="${EXON_REJOIN_MAPPING}" \
        annotation_list="${LIST_OF_ANNOTATIONS}" \
        exon_bitmasks="${EXON_BITMASKS}" \
        exon_bitmasks_coords="${EXON_BITMASK_COORDS}" \
        num_exons=${num_exons} \
        2>&1 | tee recount-unify.output.sums.txt

#now do junctions
snakemake \
    --snakefile /recount-unify/Snakefile.study_jxs \
    ${CONFIGFILE} \
    -j "${RECOUNT_CPUS}" \
    --stats recount-unify.jxs.stats.json \
    -p \
    --config \
        input=links \
        staging=staging_jxs \
        compilation_id="${PROJECT_ID}" \
        sample_ids_file="${SAMPLE_ID_MANIFEST}" \
        ref_sizes="${REF_SIZES}" \
        ref_fasta="${REF_FASTA}" \
        annotated_sjs="${ANNOTATED_JXS}" \
        build_sqlitedb=1 \
        2>&1 | tee recount-unify.output.jxs.txt
        
popd
echo SUCCESS
