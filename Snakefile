#start
import sys
import os
import glob

#current (4/13/2019) versions assumes this file structure in input:
#config[input]/study_loworder/study/run_loworder/run/groupname_attempt#
#e.g. (for CCLE, replace UUID with SRR accession if SRA/GTEx):
#ccle/le/ccle/b7/dc564d9f-3732-48ee-86ab-e21facb622b7/ccle1_in13_att2

if 'bash_path' not in config:
	bash_path='/bin/bash'

shell.executable(bash_path)

#just rejoin and after for sums
#main production version
FILES=['all.exon_bw_count.pasted.gz', 'unique.exon_bw_count.pasted.gz', 'all.logs.tar.gz', 'all.gene_counts.rejoined.tsv.gz', 'all.intron_counts.rejoined.tsv.gz', 'all.exon_counts.rejoined.tsv.gz', 'intron_counts_summed.tsv']

main_script_path=os.path.join(workflow.basedir,'scripts')

SCRIPTS={'find_done':os.path.join(main_script_path,'find_done.sh'),'find':os.path.join(main_script_path,'find_new_files.sh'),'decompress':os.path.join(main_script_path,'decompress_sums.sh'),'paste':os.path.join(main_script_path,'paste_sums.sh'),'rejoin':os.path.join(workflow.basedir, 'rejoin', 'rejoin'),'sum_counts':os.path.join(workflow.basedir, 'merge', 'sum_counts'),'QC':"python3 %s" % os.path.join(workflow.basedir, 'log_qc', 'parse_logs_for_qc.py'), 'perbase':os.path.join(workflow.basedir, 'merge', 'perbase'),'rejoin_genes':"pypy %s" % os.path.join(workflow.basedir, 'rejoin', 'rejoin_genes.py'),'split_genes':os.path.join(main_script_path,'split_out_gene_sums_by_study.sh'),'split_exons':os.path.join(workflow.basedir, 'rejoin', 'split_out_exon_sums_by_study.sh')}

#typical values for the following required parameters:
#gene_rejoin_mapping=G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed
#exon_rejoin_mapping=G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed
#gene_mapping_final=G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed
#sample_ids_file=ids.tsv
#exon_bitmasks=srav3h.exon_counts.bitmasks.tsv
#exon_bitmasks_coords=srav3h.exon_counts.bitmasks.coords
if 'gene_rejoin_mapping' not in config or 'exon_rejoin_mapping' not in config or 'sample_ids_file' not in config or 'num_samples' not in config:
	sys.stderr.write("need to pass values for 'gene_rejoin_mapping' and 'exon_rejoin_mapping' and 'sample_ids_file' and 'num_samples' for the initial rejoining part of the pipeline!\n")
	sys.exit(-1)

if 'compilation' not in config:
	sys.stderr.write("doing single study run\n")
	#sys.stderr.write("need to pass in a compilation (e.g. \"sra\" for either human or mouse, \"gtex\", or \"tcga\")\n")
	#sys.exit(-1)

#exons.bed.w_header.gz
if 'existing_sums' not in config:
	config['existing_sums']=""

gene_annotations = ["DONT_USE"]
gene_annotations_uncompressed = ["DONT_USE"]
main_annotation = None
studies = []
gene_sum_per_study_files = []
exon_sum_per_study_files = []
#order is critical, it *has* to be: "G026,G029,R109,F006,ERCC,SIRV"
#for exon splits to work correctly (bitmasks file is a static ordering of the annotations)
if 'annotation_list' in config:
	annotations_list = config['annotation_list'].split(',')
	gene_annotations_uncompressed=['%s.gene_sums.tsv' % (annotation) for annotation in annotations_list]
	gene_annotations_uncompressed.append('all.exon_counts.rejoined.tsv.gz.accession_header')
	if 'gene_mapping_final' not in config:
		sys.stderr.write("need to pass values for 'gene_mapping_final' for the final (gene) part of rejoining pipeline!\n")
		sys.exit(-1)
	gene_annotations=['%s.gene_sums.tsv.gz' % (annotation) for annotation in annotations_list]
	FILES.extend(gene_annotations)
	#typically G026
	main_annotation=annotations_list[0]
	#create FILES targets for all per-study gene & exon sums
	studies = [f.split('/')[-2:] for f in glob.glob(config['input']+'/??/*')]
	#signals we're doing a multi-study run
	if 'compilation' in config:
		if 'exon_bitmasks' not in config or 'exon_bitmasks_coords' not in config or 'num_exons' not in config:
			sys.stderr.write("need to pass values for 'num_exons' and exon_bitmasks and exon_bitmasks_coords for the final, per-study rejoining part of the pipeline!\n")
			sys.exit(-1)
		#e.g gene_sums_per_study/99/SRP214699/sra.gene_sums.SRP214699.G026.gz
		FILES.extend(["gene_sums_per_study/%s/%s/%s.gene_sums.%s.%s.gz" % (study[0], study[1], config['compilation'], study[1], annotation) for study in studies for annotation in annotations_list])
		gene_sum_per_study_files = ["gene_sums_per_study/%s/%s/%s.gene_sums.%s.{annotation}.gz" % (study[0], study[1], config['compilation'], study[1]) for study in studies]
		FILES.extend(["exon_sums_per_study/%s/%s/%s.exon_sums.%s.%s.gz" % (study[0], study[1], config['compilation'], study[1], annotation) for study in studies for annotation in annotations_list])
		exon_sum_per_study_files = ["exon_sums_per_study/%s/%s/%s.exon_sums.%s.%s.gz" % (study[0], study[1], config['compilation'], study[1], annotation) for study in studies for annotation in annotations_list]
else:
	config['gene_mapping_final'] = None
	config['annotation_list'] = None


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
		pigz --stdout -p 1 -d {input[0]} | {params.script_path} {params.gene_mapping_final_file} gene all.exon_bw_count.pasted.gz.samples_header {params.annotation_list} {params.id_mapping} {params.main_annotation}
		set +o pipefail
		paste <(echo "gene_id	chromosome	start	end	bp_length	strand") <(head -1 {params.main_annotation}.gene_sums.tsv | cut -f 6-) > all.exon_counts.rejoined.tsv.gz.accession_header
		"""

rule split_final_rejoined_genes:
	input:
		'{annotation}.gene_sums.tsv'
	output:
		gene_sum_per_study_files
	threads: 40
	params:
		script_path=SCRIPTS['split_genes'],
		compilation=config.get('compilation',''),
		annotation=lambda wildcards: wildcards.annotation
	shell:
		"""
		/bin/bash {params.script_path} {params.compilation} {params.annotation} {threads}
		"""

rule compress_final_rejoined_genes:
	input:
		'{annotation}.gene_sums.tsv'
	output:
		'{annotation}.gene_sums.tsv.gz'
	threads: 8
	params:
		annotation=lambda wildcards: wildcards.annotation
	shell:
		"""
		cat {params.annotation}.gene_sums.tsv | pigz --fast -p {threads} > {params.annotation}.gene_sums.tsv.gz
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
		cut -f 1-6 exon.counts > {output[0]}.coords
		cat exon.counts | pigz --fast -p {threads} > {output[0]}
		rm -f exon.counts exon.intron_counts
		"""

rule split_final_rejoined_exons:
	input:
		'all.exon_counts.rejoined.tsv.gz',
		'all.exon_counts.rejoined.tsv.gz.accession_header'
	output:
		exon_sum_per_study_files
	threads: 40
	params:
		script_path=SCRIPTS['split_exons'],
		compilation=config.get('compilation',''),
		annotations=config['annotation_list'],
		num_exons=config.get('num_exons',''),
		bitmasks_file=config.get('exon_bitmasks',''),
		bitmasks_coords_file=config.get('exon_bitmasks_coords','')
	shell:
		"""
		/bin/bash {params.script_path} {params.compilation} {params.annotations} {params.num_exons} {params.bitmasks_file} {params.bitmasks_coords_file} {input[0]} {threads}
        rm -rf exons_split_by_study_temp exon_annotation_split_runs
		"""

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
		config['input'],
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
