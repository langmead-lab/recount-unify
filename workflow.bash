#!/usr/bin/env bash

# Author: Chris Wilks & Ben Langmead
#  Email: broadsword@gmail.com
#   Date: 6/16/2020

# Designed to run the aggregator step (recount-unify)
# in the Monorail workflow.  This aggregates across
# multiple samples in a study, project, or compilation

set -exo pipefail

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

#TODO: use this for smoke tests
#e.g. /path/to/gene_exon_row_counts_per_annotation.tsv
#echo "Gene/exon annotation row counts file: ${GENE_EXON_ANNOTATION_ROW_COUNTS}"
#test -n "${GENE_EXON_ANNOTATION_ROW_COUNTS}"
#test -s "${GENE_EXON_ANNOTATION_ROW_COUNTS}"

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


#default case is single study, with "<sample_id>_att" style sample output directories from recount-pump
compilation_arg=""
if [[ -z $MULTI_STUDY ]]; then
    echo "Running single-study mode"
    /bin/bash -x /recount-unify/scripts/create_directory_hierarchy_for_one_study.sh $SAMPLE_ID_MANIFEST $INPUT_DIR ${WORKING_DIR}/intermediate_links > setup_intermediate_links.run 2>&1
    #need to create the directory structure we expect from the basic output of
    #recount-pump on one study (assumed)
    /bin/bash -x /recount-unify/scripts/find_done.sh ${WORKING_DIR}/intermediate_links links "*_att" > setup_links.run 2>&1
else
    echo "Running multi-study mode"
    #multi study, this assumes the input directory hierarchy/naming is correctly setup for running of find_done.sh external to recount-unify, e.g. study_loworder/study/run_loworder/run/run_inputID_att\d+
    /bin/bash -x /recount-unify/scripts/find_done.sh $INPUT_DIR links "*_att" > setup_links.run 2>&1
    compilation_arg="compilation=$PROJECT_SHORT_NAME"
fi

#e.g. 40
echo "CPUs: ${RECOUNT_CPUS}"
test -n "${RECOUNT_CPUS}"

num_samples=$(cat ${SAMPLE_ID_MANIFEST} | wc -l)
num_exons=$(cat ${EXON_BITMASK_COORDS} | wc -l)

#a config.json file could be provided in running directory
CONFIGFILE=""
if [[ -f "config.json" ]] ; then
    CONFIGFILE="--configfile config.json"
fi

#first do exon/gene sums Snakemake not skipping
if [[ -z $SKIP_SUMS ]]; then
    echo "Unifying gene and exon sums"
    snakemake \
    --snakefile /recount-unify/Snakefile \
    ${CONFIGFILE} \
    -j "${RECOUNT_CPUS}" \
    --stats recount-unify.sums.stats.json \
    -p \
    --config \
        input=links \
        staging=staging \
        sample_ids_file="${SAMPLE_ID_MANIFEST}" \
        num_samples=${num_samples} \
        existing_sums="${EXON_COORDINATES_BED}" \
        gene_rejoin_mapping="${GENE_REJOIN_MAPPING}" \
        gene_mapping_final="${GENE_ANNOTATION_MAPPING}" \
        exon_rejoin_mapping="${EXON_REJOIN_MAPPING}" \
        annotation_list="${LIST_OF_ANNOTATIONS}" \
        exon_bitmasks="${EXON_BITMASKS}" \
        exon_bitmasks_coords="${EXON_BITMASK_COORDS}" \
        num_exons=${num_exons} ${compilation_arg} \
        2>&1 | tee recount-unify.output.sums.txt

    done=`fgrep 'steps (100%) done' recount-unify.output.sums.txt`
    if [[ -z $done ]]; then
        echo "FAILURE running gene/exon unify"
        popd
        exit -1
    fi
    
    #need to do smoke tests for proper outputs based on:
    #get number of samples per-study
    #cut -f 1 $SAMPLE_ID_MANIFEST | sort | uniq -c | tr -s " " \\t | cut -f 2,3 > ${SAMPLE_ID_MANIFEST}.num_samples_per_study.tsv

    #1) check counts per study file for genes/exons
    #find gene_sums_per_study -name "*.gz" -size +0c | perl /recount-unify/scripts/check_unifier_outputs.pl ${SAMPLE_ID_MANIFEST}.num_samples_per_study.tsv genes ${GENE_EXON_ANNOTATION_ROW_COUNTS}

    #2) number of files * number of studies
    num_expected=12
    num_gene_files=$(find gene_sums_per_study -name "*.gz" -size +0c | wc -l)
    if [[ $num_expected -ne $num_gene_files ]]; then
        echo "FAILURE running gene/exon unify, unexpected # of gene sum files: $num_gene_files vs. $num_expected (expected)"
        popd
        exit -1
    fi
    num_exon_files=$(find exon_sums_per_study -name "*.gz" -size +0c | wc -l)
    if [[ $num_expected -ne $num_exon_files ]]; then
        echo "FAILURE running gene/exon unify, unexpected # of exon sum files: $num_exon_files vs. $num_expected (expected)"
        popd
        exit -1
    fi

    echo "Running QC stats collection"
    python3 /recount-unify/log_qc/parse_logs_for_qc.py --incoming-dir links --sample-mapping "${SAMPLE_ID_MANIFEST}" --intron-sums intron_counts_summed.tsv > qc_1.tsv 2> qc.err
    num_samples_qc=$(cat qc_1.tsv | wc -l)
    num_samples_qc=$(( num_samples_qc - 1 ))
    if [[ $num_samples_qc -ne $num_samples ]]; then
        echo "FAILURE in pump output QC stats collection (post gene/exon unify), # QC rows ($num_samples_qc) != # samples ($num_samples)"
        popd
        exit -1
    fi
fi

#now do junctions if not skipping
if [[ -z $SKIP_JUNCTIONS ]]; then
    echo "Unifying junction counts"
    #need to override the default python path (conda python3)
    #to get python2.7 with PyLucene installed, but only for build_lucene purposes
    export PYPATH=/usr/bin/python
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
        sample_original_metadata_file="${SAMPLE_ID_MANIFEST_ORIG}" \
        ref_sizes="${REF_SIZES}" \
        ref_fasta="${REF_FASTA}" \
        annotated_sjs="${ANNOTATED_JXS}" \
        study_dir="junction_counts_per_study" \
        build_sqlitedb=1 \
        build_lucene=1 ${compilation_arg} \
        2>&1 | tee recount-unify.output.jxs.txt

    done=`fgrep 'steps (100%) done' recount-unify.output.jxs.txt`
    if [[ -z $done ]]; then
        echo "FAILURE running junction unify"
        popd
        exit -1
    fi
    #do a little clean up not part of recount-unify proper
    if [[ ! -z $MULTI_STUDY ]]; then
        mkdir -p temp_jxs
        mv junction_counts_per_study/?? temp_jxs/
        mv junction_counts_per_study junction_counts_per_study_run_files
        mv temp_jxs junction_counts_per_study
    fi
    num_jx_files=$(find junction_counts_per_study -name "*.gz" -size +0c | wc -l)
    if [[ $num_expected -ne $num_jx_files ]]; then
        echo "FAILURE running gene/exon unify, unexpected # of gene sum files: $num_jx_files vs. $num_expected (expected)"
        popd
        exit -1
    fi

    #add jx stats to QC stats
    cat qc_1.tsv | perl /recount-unify/log_qc/add_jx_stats2qc.pl samples.tsv > qc_2.tsv 2>> qc.err
    #this step adds in the combined columns for the STAR stats in the case of 2-step alignment when there's a 3rd, unpaired FASTQ file
    python3 /recount-unify/log_qc/add_together_qc_fields.py qc_2.tsv > qc.tsv 2>> qc.err
    #do a 2nd check here at the end of the QC stats
    num_samples_qc=$(cat qc.tsv | wc -l)
    num_samples_qc=$(( num_samples_qc - 1 ))
    if [[ $num_samples_qc -ne $num_samples ]]; then
        echo "FAILURE in pump output QC stats collection (post gene/exon unify), # QC rows ($num_samples_qc) != # samples ($num_samples)"
        popd
        exit -1
    fi
    popd
fi

echo SUCCESS
