#start
import sys
import os
import glob

#current (4/13/2019) versions assumes this file structure in input:
#config[input]/study_loworder/study/run_loworder/run/groupname_attempt#
#e.g. (for CCLE, replace UUID with SRR accession if SRA/GTEx):
#ccle/le/ccle/b7/dc564d9f-3732-48ee-86ab-e21facb622b7/ccle1_in13_att2

#just rejoin and after for sums
#FILES=['all.gene_counts.rejoined.tsv.gz', 'all.intron_counts.rejoined.tsv.gz', 'all.exon_counts.rejoined.tsv.gz', 'intron_counts_summed.tsv']
#main production version
FILES=['all.exon_bw_count.pasted.gz', 'unique.exon_bw_count.pasted.gz', 'all.sjs.merged.annotated.tsv.gz', 'all.logs.tar.gz', 'all.gene_counts.rejoined.tsv.gz', 'all.intron_counts.rejoined.tsv.gz', 'all.exon_counts.rejoined.tsv.gz', 'intron_counts_summed.tsv']

main_script_path=os.path.join(workflow.basedir,'scripts')

SCRIPTS={'find_done':os.path.join(main_script_path,'find_done.sh'),'find':os.path.join(main_script_path,'find_new_files.sh'),'decompress':os.path.join(main_script_path,'decompress_sums.sh'),'paste':os.path.join(main_script_path,'paste_sums.sh'),'filter':os.path.join(main_script_path,'filter_new_sjs.sh'),'merge':os.path.join(workflow.basedir, 'merge', 'merge.py'),'annotate':os.path.join(workflow.basedir, 'annotate', 'annotate_sjs.py'),'rejoin':os.path.join(workflow.basedir, 'rejoin', 'rejoin'),'sum_counts':os.path.join(workflow.basedir, 'merge', 'sum_counts'),'QC':"python3 %s" % os.path.join(workflow.basedir, 'log_qc', 'parse_logs_for_qc.py'), 'perbase':os.path.join(workflow.basedir, 'merge', 'perbase'),'rejoin_genes':"pypy %s" % os.path.join(workflow.basedir, 'rejoin', 'rejoin_genes.py')}

#typical values for the following required parameters:
#gene_rejoin_mapping=G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed
#exon_rejoin_mapping=G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed
#gene_mapping_final=G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed
#sample_ids_file=ids.tsv
if 'gene_rejoin_mapping' not in config or 'exon_rejoin_mapping' not in config or 'gene_mapping_final' not in config or 'sample_ids_file' not in config or 'num_samples' not in config:
	sys.stderr.write("need to pass values for 'gene_rejoin_mapping' and/or 'exon_rejoin_mapping' and/or 'gene_mapping_final' and/or 'sample_ids_file' and/or 'num_samples' for the rejoining part of the pipeline!\n")
	sys.exit(-1)	

#ref_sizes=hg38.recount_pump.fa.new_sizes
#ref_fasta=hg38.recount_pump.fa
if 'ref_sizes' not in config or 'ref_fasta' not in config:
	sys.stderr.write("need to pass values for 'ref_sizes' and/or 'ref_fasta' the jx motif extraction part of the pipeline!\n")
	sys.exit(-1)	

#annotated_junctions.tsv.gz
if 'existing_sj_db' not in config:
	config['existing_sj_db']=""

#exons.bed.w_header.gz
if 'existing_sums' not in config:
	config['existing_sums']=""

#e.g. G026,G029,R109,ERCC,SIRV,F006
if 'annotation_list' not in config:
	sys.stderr.write("need to pass value(s) for 'annotation_list' (e.g. 'G026,G029,R109,ERCC,SIRV,F006') for final gene sums!\n")
	sys.exit(-1)

gene_annotations=['%s.gene.sums.tsv.gz' % (annotation) for annotation in config['annotation_list'].split(',')]
FILES.extend(gene_annotations)
#typically G026
main_annotation=config['annotation_list'].split(',')[0]

if 'compilation_id' not in config:
	config['compilation_id']=0

wildcard_constraints:
	study_group_num="[0-9a-zA-Z]{2}",
	run_group_num="[0-9a-zA-Z]{2}",
	type="(all)|(unique)"

rule all:
	input:
		expand("{file}", file=FILES)

#this *HAS* to be outside of and before Snakemake, since
#the wildcard groupings can't be inferred w/o it being run first
#otherwise the rule dependency graph will skip critical steps (e.g. find_sums)
#rule find_done:
#	input:
#		config['recount_pump_output']
#	output:
#		config['input']
#	params:
#		input_dir=config['input'],
#		script_path=SCRIPTS['find_done'],
#		study=config['study']
##	shell:
#		"""
#		{params.script_path} {input} {params.input_dir} {params.study}
#		""" 


#tar and gzip all the logs, but maintain the directory structure
rule tar_logs:
	input:
		config['input']
	output:
		'all.logs.tar.gz'
	shell:
		"find -L {input} -name '*.log' > all_logs && tar -zcvf {output} -T all_logs > /dev/null && rm all_logs"


###exon SUM pasting rules
rule find_sums:
	input: 
		config['input'],
		config['sample_ids_file']
	output:
		config['staging'] + '/{type}.exon_bw_count.groups.manifest'
	params:
		staging=config['staging'],
		script_path=SCRIPTS['find'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {input[0]} {input[1]} {params.staging} {params.type}.exon_bw_count .zst"

rule decompress_sums:
	input:
		config['staging'] + '/{type}.exon_bw_count.groups.manifest'
	output:
		config['staging'] + '/{type}.exon_bw_count.{study_group_num}.{run_group_num}.decompressed'
	params:
		study_group_num=lambda wildcards: wildcards.study_group_num,
		run_group_num=lambda wildcards: wildcards.run_group_num,
		staging=config['staging'],
		script_path=SCRIPTS['decompress'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {params.staging}/{params.type}.exon_bw_count.{params.study_group_num}.{params.run_group_num}.manifest {output}"

#do a rule instantiation per *run* low-order name grouping to do hierarchical pastes
rule paste_sums_per_group:
	input:
		config['staging'] + '/{type}.exon_bw_count.{study_group_num}.{run_group_num}.decompressed'
	output:
		config['staging'] + '/{type}.exon_bw_count.{study_group_num}.{run_group_num}.pasted'
	params:
		study_group_num=lambda wildcards: wildcards.study_group_num,
		run_group_num=lambda wildcards: wildcards.run_group_num,
		staging=config['staging'],
		script_path=SCRIPTS['paste'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {params.staging}/{params.type}.exon_bw_count.{params.study_group_num}.{params.run_group_num}.manifest {output}"

def get_pasted_sum_files(wildcards):
	study_loworder = wildcards.study_group_num
	return [config['staging']+"/%s.exon_bw_count.%s.%s.pasted" % (wildcards.type, f.split('/')[-3], f.split('/')[-1]) for f in glob.glob(config['input']+'/%s/*/??' % (study_loworder))]

rule collect_pasted_sums:
	input:
		get_pasted_sum_files
	output:
		config['staging'] + '/{type}.exon_bw_count.{study_group_num}.pasted.files.list'
	params:
		study_group_num=lambda wildcards: wildcards.study_group_num,
		staging=config['staging'],
		type=lambda wildcards: wildcards.type
	shell:
		"ls {params.staging}/{params.type}.exon_bw_count.{params.study_group_num}.??.pasted > {output}"

rule paste_sums_per_study_group:
	input:
		config['staging'] + '/{type}.exon_bw_count.{study_group_num}.pasted.files.list'
	output:
		os.path.join(config['staging'], '{type}.exon_bw_count.{study_group_num}.pasted')
	params:
		study_group_num=lambda wildcards: wildcards.study_group_num,
		staging=config['staging'],
		script_path=SCRIPTS['paste'],
		existing_sums=config['existing_sums'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {input} {output} dont_get_ids"


def get_study_pasted_sum_files(wildcards):
	return [config['staging']+"/%s.exon_bw_count.%s.pasted" % (wildcards.type, f.split('/')[-1]) for f in glob.glob(config['input']+'/??')]

rule collect_study_pasted_sums:
	input:
		get_study_pasted_sum_files
	output:
		config['staging'] + '/{type}.exon_bw_count.groups.pasted.files.list'
	params:
		staging=config['staging'],
		type=lambda wildcards: wildcards.type
	shell:
		"ls {params.staging}/{params.type}.exon_bw_count.??.pasted > {output}"

rule paste_sums_final:
	input:
		config['staging'] + '/{type}.exon_bw_count.groups.pasted.files.list'
	output:
		'{type}.exon_bw_count.pasted.gz'
	params:
		staging=config['staging'],
		script_path=SCRIPTS['paste'],
		existing_sums=config['existing_sums'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {input} {output} dont_get_ids {params.existing_sums}"


###Rejoin of exon/gene coverage into original annotation rows
rule rejoin_genes:
	input:
		'all.exon_bw_count.pasted.gz'
	output:
		'all.gene_counts.rejoined.tsv.gz',
		'all.intron_counts.rejoined.tsv.gz'
	threads: 8
	params:
		staging=config['staging'],
		script_path=SCRIPTS['rejoin'],
		gene_mapping_file=config['gene_rejoin_mapping'],
		num_samples=config['num_samples']
	shell:
		"""
		{params.script_path} -a {params.gene_mapping_file} -d <(pigz --stdout -p 1 -d {input}) -s {params.num_samples} -p gene -h  
	        cat gene.counts | pigz --fast -p {threads} > {output[0]}
		cat gene.intron_counts | pigz --fast -p {threads} > {output[1]}
		rm -f gene.counts gene.intron_counts
		"""

gene_annotations_uncompressed=['%s.gene.sums.tsv' % (annotation) for annotation in config['annotation_list'].split(',')]
rule rejoin_genes_final:
	input:
		'all.gene_counts.rejoined.tsv.gz',
		'all.exon_bw_count.pasted.gz'
	output:
		gene_annotations_uncompressed
	params:
		staging=config['staging'],
		script_path=SCRIPTS['rejoin_genes'],
		#../G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed
		gene_mapping_final_file=config['gene_mapping_final'],
		#G026,G029,R109,ERCC,SIRV,F006
		annotation_list=config['annotation_list'],
		id_mapping=config['sample_ids_file'],
        	#G026
		main_annotation=main_annotation
	shell:
		"""
		set +o pipefail
		pigz --stdout -p 1 -d {input[1]} | head -1 | cut -f 7- > all.exon_bw_count.pasted.gz.samples_header
		set -o pipefail
		pigz --stdout -p 1 -d {input[0]} | tail -n+2 | {params.script_path} {params.gene_mapping_final_file} gene all.exon_bw_count.pasted.gz.samples_header {params.annotation_list} {params.id_mapping} {params.main_annotation}
		"""

rule compress_final_rejoined_genes:
	input:
		'{annotation}.gene.sums.tsv'
	output:
		'{annotation}.gene.sums.tsv.gz'
	threads: 8
	params:
		annotation=lambda wildcards: wildcards.annotation
	shell:
		"""
		pigz --fast -p {threads} {params.annotation}.gene.sums.tsv
		"""

rule rejoin_exons:
	input:
		'all.exon_bw_count.pasted.gz'
	output:
		'all.exon_counts.rejoined.tsv.gz'
	threads: 8
	params:
		staging=config['staging'],
		script_path=SCRIPTS['rejoin'],
		exon_mapping_file=config['exon_rejoin_mapping'],
		num_samples=config['num_samples']
	shell:
		"""
		set +o pipefail
		{params.script_path} -a {params.exon_mapping_file} -d <(pigz --stdout -p 1 -d {input}) -s {params.num_samples} -p exon -h 
	        cat exon.counts | pigz --fast -p {threads} > {output[0]}
		rm -f exon.counts exon.intron_counts
		"""

###Splice junction merging rules
rule find_sjs:
	input: 
		config['input'],
		config['sample_ids_file']
	output:
		config['staging'] + '/sj.groups.manifest'
	params:
		staging=config['staging'],
		wildc='"*.zst"',
		script_path=SCRIPTS['find']
	shell:
		"{params.script_path}  {input[0]} {input[1]} {params.staging} sj {params.wildc}"

rule filter_sjs:
	input:
		config['staging'] + '/sj.groups.manifest'
	output:
		config['staging'] + '/sj.{study_group_num}.{run_group_num}.manifest.filtered'
	params:
		study_group_num=lambda wildcards: wildcards.study_group_num,
		run_group_num=lambda wildcards: wildcards.run_group_num,
		staging=config['staging'],
		script_path=SCRIPTS['filter']
	shell:
		"{params.script_path} {params.staging}/sj.{params.study_group_num}.{params.run_group_num}.manifest"

#merge at the group level, uses sample_ids
rule merge_sjs:
	input:
		config['staging'] + '/sj.{study_group_num}.{run_group_num}.manifest.filtered'
	output:
		config['staging'] + '/sj.{study_group_num}.{run_group_num}.merged.tsv.gz'
	threads: 8
	params:
		staging=config['staging'],
		filtered_manifest=lambda wildcards, input: '.'.join(input[0].split('.')[:-1]),
		script_path=SCRIPTS['merge']
	shell:
		"pypy {params.script_path} --list-file {params.filtered_manifest} --gzip | sort -k1,1 -k2,2n -k3,3n | pigz --fast -p {threads} > {output}"

#gets  study + run loworder groupings, e.g. unified/sj.01.42.manifest.sj_sample_files.merged.tsv.gz
#where "42" is the loworder digits for the run, e.g. SRR8733242 in SRP005401
def get_sj_merged_files(wildcards):
	study_loworder = wildcards.study_group_num
	return [config['staging']+"/sj.%s.%s.merged.tsv.gz" % (f.split('/')[-3],f.split('/')[-1]) for f in glob.glob(config['input']+'/%s/*/??' % (study_loworder))]

rule collect_merged_sjs:
	input:
		get_sj_merged_files
	output:
		config['staging'] + '/sj.{study_group_num}.groups.merged.files.list'
	params:
		staging=config['staging'],
		study_group_num=lambda wildcards: wildcards.study_group_num
	shell:
		"ls {params.staging}/sj.{params.study_group_num}.??.merged.tsv.gz > {output}"

rule merge_study_group_sjs:
	input: 
		config['staging'] + '/sj.{study_group_num}.groups.merged.files.list'
	output:
		config['staging'] + '/all.{study_group_num}.sjs.merged.tsv.gz'
	threads: 8
	params:
		script_path=SCRIPTS['merge'],
		existing_sj_db=config['existing_sj_db']
	shell:
		"pypy {params.script_path} --list-file {input} --gzip --append-samples | pigz --fast -p {threads} > {output}"

def get_sj_study_merged_files(wildcards):
	return [config['staging']+"/all.%s.sjs.merged.tsv.gz" % (f.split('/')[-1]) for f in glob.glob(config['input']+'/??')]

rule collect_study_merged_sjs:
	input:
		get_sj_study_merged_files
	output:
		config['staging'] + '/sj.groups.merged.files.list'
	params:
		staging=config['staging']
	shell:
		"ls {params.staging}/all.??.sjs.merged.tsv.gz > {output}"

rule merge_all_sjs:
	input: 
		config['staging'] + '/sj.groups.merged.files.list'
	output:
		config['staging'] + '/all.sjs.merged.tsv'
	params:
		script_path=SCRIPTS['merge'],
		existing_sj_db=config['existing_sj_db']
	shell:
		"pypy {params.script_path} --list-file {input} --gzip --append-samples --existing-sj-db \"{params.existing_sj_db}\" > {output}"

rule extract_motifs_for_sjs:
	input:
		config['staging'] + '/all.sjs.merged.tsv',
		config['ref_sizes'],
		config['ref_fasta']
	output:
		config['staging'] + '/all.sjs.merged.motifs.tsv'
	params:
		script_path=SCRIPTS['perbase'],
	shell:
		"""
                cat {input[0]} | {params.script_path} -c {input[1]} -g {input[2]} -f {input[1]} > {output}.temp 2>{output}.errrs
		paste <(cut -f 1-4 {output}.temp) <(cut -f 7 {output}.temp | tr '[:lower:]' '[:upper:]') <(cut -f 6 {output}.temp) > {output}
		"""
#| perl -ne '$f=$_; @f=split(/\t/,$_); $f[4]=uc($f[4]); print "".join("\t",@f)."\n";' > {output}

rule annotate_all_sjs:
	input:
		config['staging'] + '/all.sjs.merged.motifs.tsv'
	output:
		'all.sjs.merged.annotated.tsv.gz'
	threads: 8
	params:
		annot_sjs=config['annotated_sjs'],
		script_path=SCRIPTS['annotate'],
		compilation_id=config['compilation_id']
	shell:
		"cat {input} | pypy {params.script_path} --compiled-annotations {params.annot_sjs} --compilation-id {params.compilation_id} | pigz --fast -p {threads} > {output}"


rule sum_intron_counts:
	input:
		'all.exon_bw_count.pasted.gz',
		'all.intron_counts.rejoined.tsv.gz'
	output:
		'intron_counts_summed.tsv'
	params:
		script_path=SCRIPTS['sum_counts']
	shell:
		"""
		set +eo pipefail
		pigz --stdout -p 1 -d {input[0]} | cut -f 6- | head -1 > {output}
		set -eo pipefail
		pigz --stdout -p 1 -d {input[1]} | {params.script_path} >> {output} 2> {output}.err
		"""

#while this can be run within snakemake
#it's more modular to keep it outside
rule QC:
	input:
		config['recount_pump_output'],
		'intron_counts_summed.tsv'
	output:
		'qc.stats.tsv'
	params:
		script_path=SCRIPTS['QC'],
		sample_ids=config['sample_ids_file']
	shell:
		"""
		{params.script_path} --incoming-dir {input[0]} --sample-mapping {params.sample_ids} --intron-sums {input[1]} > {output} 2> {output}.err
		"""
