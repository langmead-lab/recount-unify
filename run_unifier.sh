#!/usr/bin/bash
#run the whole recount unifier pipeline
#1) determine all finished (done) runs and symlink them in expected hierarchy
#2) run Snakemake summarization pipeline on output of 1)
#3) run QC stats on original file paths but requires the intron sums from 2)

#NOTE: script assumes:
#1) it's being run from a tranche/compilation specific subdirectory
#   under the recount-unifier repo checkout root directory
#2) the total number of unique run accessions/uuids successfully processed is known (2nd argument below)

#may need to update these with different paths/filenames/fileversions

##These files are for getting the Snaptron-formatted junction coverage files
#used to determine which junctions calls coming from Monorail are previously annotated and from which annotation
annotated_sjs='../annotated_junctions.tsv.gz'
#need the reference genome to determine the motifs of the splice sites (efficiently)
ref_sizes='../hg38.recount_pump.fa.new_sizes'
ref_fasta='../hg38.recount_pump.fa'

##These files are used for summarizing the sums back to their original annotation intervals (genes/exons)
#this is the coordinates of every disjoint exon we use, it *has* to be the same order
#as the annotation file used to generate the exon sums
existing_sums='../exons.bed.w_header.gz'
#these 4 files are used in the "rejoin" process to get sums
#for the original individual annotations' genes/exons intervals
gene_rejoin_mapping='../G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed'
exon_rejoin_mapping='../G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed'
gene_mapping_final='../G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed'
annotation_list='G026,G029,R109,ERCC,SIRV,F006'

#compilation ID used in rail_id generation (e.g. tcga is 3, sra_human_v3_8 is 18, sra_human_v3_2 is 12)
comp_id=$1
#number of run accessions (SRRs), this is the set of unique runs that were actually processed
#which almost always be a true subset of the input set of runs (contained in ids.tsv, below)
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

#example run for TCGA
#/bin/bash -x ../run_unifier.sh 3 11373 40 /storage2/cwilks/recount_pump/destination/tcga tcga
#or SRA human tranche 0:
#/bin/bash -x ../run_unifier.sh 10 30000 40 /storage2/cwilks/recount_pump/destination/srav3/sra_human_v3_0 sra_human_v3

python2 ../sample_ids/assign_compilation_ids.py --sep ' ' --acc-col 1 --accessions-file $accessions_file --compilation-code $comp_id > sample_rail_ids.${comp_id}

ln -fs sample_rail_ids.${comp_id} ids.tsv

ln -fs ../blank_exon_sums
ln -fs $original_full_path

/bin/bash -x ../scripts/find_done.sh $original_full_path links $study_prefix

timebin=`which time`

$timebin -v snakemake -j $threads --stats ./stats.json --snakefile ../Snakefile -p --config input=links staging=unified sample_ids_file=ids.tsv annotated_sjs=$annotated_sjs existing_sums=$existing_sums compilation_id=$comp_id gene_rejoin_mapping=$gene_rejoin_mapping exon_rejoin_mapping=$exon_rejoin_mapping num_samples=$num_runs ref_sizes=$ref_sizes ref_fasta=$ref_fasta gene_mapping_final=$gene_mapping_final annotation_list=$annotation_list > summarize.run 2>summarize.err

python3 ../log_qc/parse_logs_for_qc.py --incoming-dir links --sample-mapping ids.tsv --intron-sums intron_counts_summed.tsv > qc${comp_id}.tsv 2> qc${comp_id}.err
