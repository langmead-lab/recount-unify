#!/usr/bin/bash
#run the whole recount unifier pipeline
#1) determine all finished (done) runs and symlink them in expected hierarchy
#2) run gene/exon summarization pipeline ("Snakefile") on output of 1)
#4) run junction summarization pipeline ("Snakefile.study_jxs") on output of 1)
#5) run QC stats on original file paths but requires the intron sums from 2)

#NOTE: script assumes:
#1) it's being run from a tranche/compilation specific subdirectory
#   under the recount-unifier repo checkout root directory
#2) the total number of unique run accessions/uuids successfully processed is known (2nd argument below), this has to be exact

dir=$(dirname $0)

#compilation ID used in rail_id generation (e.g. tcga is 3, sra_human_v3_8 is 18, sra_human_v3_2 is 12)
comp_id=$1
#number of run accessions (SRRs), this is the set of unique runs that were actually processed
#which will almost always be a true subset of the input set of runs (contained in ids.tsv, below)
#it also *does not* include the 1) "spike-ins" and the 2) overlapping studies from the next tranche (those would just be dups)
num_runs=$2
threads=$3
#e.g. /storage2/cwilks/recount_pump/destination/srav3/sra_human_v3_9
original_full_path=$4
#e.g. tcga,sra_human_v3,etc...
study_prefix=$5
#this is the full path to the sample accessions in the project/tranche (used to generate the sample/rail_ids)
#this is a space delimited file with the run accessions/uuids as the 2nd column
accessions_file=$6
#either "hg38" or "grcm38"
orgn=$7
#annotation list, e.g. "G026,G029,R109,F006,ERCC,SIRV" (human) or "M023,ERCC,SIRV" (mouse)
annotation_list=$8
#path to where all the unifier genome indexes and related are, e.g. ../refs
REFS=$9

REFS="$REFS/${orgn}_unify"

##These files are for getting the Snaptron-formatted junction coverage files
#used to determine which junctions calls coming from Monorail are previously annotated and from which annotation
annotated_sjs="$REFS/annotated_junctions.tsv.gz"
#this should be the exact same FASTA file used to generate the genome reference indexes in the pump runs
ref_fasta="$REFS/recount_pump.fa"
#need the reference genome to determine the motifs of the splice sites (efficiently)
ref_sizes="$REFS/recount_pump.chr_sizes.tsv"

##These files are used for summarizing the sums back to their original annotation intervals (genes/exons)
#this is the coordinates of every disjoint exon we use, it *has* to be the same order
#as the annotation file used to generate the exon sums
existing_sums="$REFS/exons.w_header.bed.gz"
num_exons=$(pcat $existing_sums | tail -n+2 | wc -l)

exon_bitmasks="$REFS/exon_bitmasks.tsv"
exon_bitmasks_coords="$REFS/exon_bitmask_coords.tsv"

#these 4 files are used in the "rejoin" process to get sums
#for the original individual annotations' genes/exons intervals
gene_rejoin_mapping="$REFS/disjoint2exons2genes.bed"
exon_rejoin_mapping="$REFS/disjoint2exons.bed"
gene_mapping_final="$REFS/disjoint2exons2genes.rejoin_genes.bed"

#example run for TCGA
#/bin/bash -x ../run_unifier.sh 3 11373 40 /storage2/cwilks/recount_pump/destination/tcga tcga
#or SRA human tranche 0:
#/bin/bash -x ../run_unifier.sh 10 30000 40 /storage2/cwilks/recount_pump/destination/srav3/sra_human_v3_0 sra_human_v3

python2 $dir/sample_ids/assign_compilation_ids.py --sep ' ' --acc-col 1 --accessions-file $accessions_file --compilation-code $comp_id > sample_rail_ids.${comp_id}

ln -fs sample_rail_ids.${comp_id} ids.tsv

ln -fs $REFS/blank_exon_sums ./blank_exon_sums
ln -fs $original_full_path

/bin/bash -x $dir/scripts/find_done.sh $original_full_path links $study_prefix

timebin=`which time`

$timebin -v snakemake -j $threads --stats ./stats.json --snakefile $dir/Snakefile -p --config input=links staging=unified sample_ids_file=ids.tsv annotated_sjs=$annotated_sjs existing_sums=$existing_sums compilation=$comp gene_rejoin_mapping=$gene_rejoin_mapping exon_rejoin_mapping=$exon_rejoin_mapping gene_mapping_final=$gene_mapping_final num_samples=$num_runs num_exons=$num_exons annotation_list=$annotation_list > gene_exon_summarize.run 2> gene_exon_summarize.err

#compilation_id=$comp_id ref_sizes=$ref_sizes ref_fasta=$ref_fasta

python3 $dir/log_qc/parse_logs_for_qc.py --incoming-dir links --sample-mapping ids.tsv --intron-sums intron_counts_summed.tsv > qc${comp_id}.tsv 2> qc${comp_id}.err
