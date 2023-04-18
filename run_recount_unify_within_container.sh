#!/usr/bin/env bash
set -exo pipefail

#For dbGaP runs, define PROTECTED=1 in the running environment (e.g. in bash "export PROTECTED=1")
#To do genotyping with QC, define SNPS_FILE_FOR_GENOTYPING=$REF_DIR_HOST/path/to/snps.vcf
#To bypass initial link creation for a rerun, define SKIP_FIND=1

#hg38 or grcm38
export REF=$1

#full path to one directory above where get_unify_refs.sh deposited its files
export REF_DIR=$2

#full path to where we should actually run the unifier
export WORKING_DIR=$3

#full path to where the output from recount-pump resides
#this needs to be writable by this script!
export INPUT_DIR_HOST=$4

#list of 2 or more study_id<TAB>sample_id + any number of tab delimited optional fields
#this file MUST have a header otherwise the first row will be skipped!
export SAMPLE_ID_MANIFEST_HOST=$5

#number of processes to start within container, 10-40 are reasonable depending on the system/run
export NUM_CPUS=$6

#optional, this is used as the project name as well as compilation_id in the jx output, defaults to "sra" and 101 respectively
#compilation short name (e.g. "sra", "gtex", or "tcga", to be compatible with recount3 outputs) or custom name
#format: 'project_short_name:project_id', default:'sra:101'
export PROJECT_SHORT_NAME='sra'
export PROJECT_ID=101
export PROJECT_SHORT_NAME_AND_ID=$8
if [[ -n $PROJECT_SHORT_NAME_AND_ID ]]; then
    failed_format_check=$(perl -e '$in="'$PROJECT_SHORT_NAME_AND_ID'"; chomp($in); ($p,$pid)=split(/:/,$in); if(length($p) == 0 || $pid!~/^(\d+)$/ || $pid < 100 || $pid > 249) {  print" bad project short name ($p) and/or project ID ($pid) input, format <project_short_name(str)>:project_id(int)> and project_id must be > 100 and < 250, exiting\n"; exit(-1);}')
    if [[ -n $failed_format_check ]]; then
        echo $failed_format_check
        exit -1
    fi
    export PROJECT_SHORT_NAME=$(echo $PROJECT_SHORT_NAME_AND_ID | tr ':' \\t | cut -f 1) 
    export PROJECT_ID=$(echo $PROJECT_SHORT_NAME_AND_ID | tr ':' \\t | cut -f 2) 
fi
echo "PROJECT_SHORT_NAME=$PROJECT_SHORT_NAME"
echo "PROJECT_ID=$PROJECT_ID"

export MULTI_STUDY=1
#optional, only used if you explicitly want recount-unify to build a single study
#this is usually true only when you want to skip producing recount3 formatted data
#e.g. you only want Snaptron-ready output
SINGLE_STUDY_ONLY=$9
if [[ ! -z $SINGLE_STUDY_ONLY ]]; then
    MULTI_STUDY=
fi

export INPUT_DIR=$WORKING_DIR/input_from_pump
mkdir -p $INPUT_DIR

#NOTE: the following assumes the unifier is being run on the *same filesystem* as the pump
#as it makes hard links to the pump output
#make sure input data is properly organized for the Unifier
#assumes an original output format: $INPUT_DIR_HOST/sampleID_att0/sampleID!studyID!*.manifest
#we can skip this if $SKIP_FIND is set in the running environment
#../geuvadis_small_output/ERR188431_att0/ERR188431!ERP001942!hg38!sra.align.log
if [[ ! -z $MULTI_STUDY  && -z $SKIP_FIND ]]; then
    find $INPUT_DIR_HOST -name '*!*' | perl -ne 'BEGIN { $run_id=1; } $working_dir="'$INPUT_DIR'"; chomp; $f=$_; @f=split(/\//,$f); $fm=pop(@f); $original=join("/",@f); $run_dir=pop(@f); @f2=split(/!/,$fm); $sample=shift(@f2); if(!$h{$sample}) { $h{$sample}=$run_id++; } $i=$h{$sample}; $study=shift(@f2); $study=~/(..)$/; $lo1=$1; $sample=~/(..)$/; $lo2=$1; $parent=join("/",@f); $newsub="$working_dir/$lo1/$study/$lo2/$sample"; $i++; $run_dir=~s/(_att\d+)$/_in$i$1/;  `mkdir -p $newsub/$run_dir ; ln -f $f $newsub/$run_dir/ ; touch $newsub/$run_dir.done`;'
fi

#TODO: work out *within container* directories in this new version
export REF_DIR=$REF_DIR/$REF"_unify"

#human
if [[ $REF == 'hg38' ]]; then
    export LIST_OF_ANNOTATIONS='G026,G029,R109,F006,ERCC,SIRV'
#mouse
else 
    if [[ $REF == 'grcm38' ]]; then
        export LIST_OF_ANNOTATIONS='M023,ERCC,SIRV'
    else
        echo "ERROR: unrecognized organism: $REF, exiting"
        exit
    fi
fi

#generic names are used for the organism specific REF files upstream
#so just need to assign them to the ENV vars expected by the container
export ANNOTATED_JXS=$REF_DIR/annotated_junctions.tsv.gz
export EXON_COORDINATES_BED=$REF_DIR/exons.w_header.bed.gz
export EXON_REJOIN_MAPPING=$REF_DIR/disjoint2exons.bed
export GENE_REJOIN_MAPPING=$REF_DIR/disjoint2exons2genes.bed
export GENE_ANNOTATION_MAPPING=$REF_DIR/disjoint2exons2genes.rejoin_genes.bed
export REF_FASTA=$REF_DIR/recount_pump.fa
export REF_SIZES=$REF_DIR/recount_pump.chr_sizes.tsv
export EXON_BITMASKS=$REF_DIR/exon_bitmasks.tsv
export EXON_BITMASK_COORDS=$REF_DIR/exon_bitmask_coords.tsv

#export INPUT_DIR=$INPUT_DIR
#export WORKING_DIR=$WORKING_DIR

export RECOUNT_CPUS=$NUM_CPUS


#do some munging of the passed in sample IDs + optional metadata files
sample_id_manfest_fn=$(basename $SAMPLE_ID_MANIFEST_HOST)
#first copy the full original samples manifest into a directory visible to the container
cp $SAMPLE_ID_MANIFEST_HOST $WORKING_DIR/$sample_id_manfest_fn
export SAMPLE_ID_MANIFEST_HOST_ORIG=$WORKING_DIR/$sample_id_manfest_fn
export SAMPLE_ID_MANIFEST_HOST=$WORKING_DIR/ids.input
#now cut out just the 1st 2 columns (study, sample_id), expecting a header
cut -f 1,2 $SAMPLE_ID_MANIFEST_HOST_ORIG | tail -n+2 > $SAMPLE_ID_MANIFEST_HOST
export SAMPLE_ID_MANIFEST=$WORKING_DIR/ids.input
export SAMPLE_ID_MANIFEST_ORIG=$WORKING_DIR/$sample_id_manfest_fn

export ORGANISM_REF=$REF

#perl -e 'print "SAMPLE_ID_MANIFEST\nREF_DIR\nSAMPLE_ID_MANIFEST_ORIG\nRECOUNT_CPUS\nWORKING_DIR\nINPUT_DIR\nEXON_BITMASK_COORDS\nEXON_BITMASKS\nREF_SIZES\nREF_FASTA\nGENE_ANNOTATION_MAPPING\nGENE_REJOIN_MAPPING\nEXON_REJOIN_MAPPING\nEXON_COORDINATES_BED\nANNOTATED_JXS\nLIST_OF_ANNOTATIONS\nPROJECT_SHORT_NAME\nPROJECT_ID\nORGANISM_REF\nMULTI_STUDY\nSNPS_FILE_FOR_GENOTYPING\nPROTECTED";' > $WORKING_DIR/env.list

#docker run --rm --user $RU_UID --name runifer --env-file $WORKING_DIR/env.list --volume $INPUT_DIR:$INPUT_DIR --volume $WORKING_DIR:${WORKING_DIR}:rw --volume $REF_DIR_HOST:${REF_DIR}:rw $image_name /bin/bash -c "/recount-unify/workflow.bash.recount_jxns_only"
/usr/bin/time -v /bin/bash -x /recount-unify/workflow.bash.recount_jxns_only

#putting all relevant final output files in one directory
mkdir -p $WORKING_DIR/run_files
pushd $WORKING_DIR
ls | fgrep -v run_files | xargs -I {} mv {} run_files/
pushd $WORKING_DIR/run_files
mv ids.tsv samples.* gene_sums_per_study exon_sums_per_study junction_counts_per_study metadata $WORKING_DIR/
popd
popd
