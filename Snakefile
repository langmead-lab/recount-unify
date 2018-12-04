#start
import os
import glob

FILES=[os.path.join(config['staging'], config['study'] + '.all.sjs.merged.annotated.tsv.gz'), os.path.join(config['staging'], config['study'] + '.all.exon_bw_count.pasted.gz'), os.path.join(config['staging'], config['study'] + '.unique.exon_bw_count.pasted.gz'), os.path.join(config['staging'], config['study'] + '.all.logs.tar.gz')]
main_script_path=os.path.join(workflow.basedir,'scripts')
SCRIPTS={'find':os.path.join(main_script_path,'find_new_files.sh'),'decompress':os.path.join(main_script_path,'decompress_sums.sh'),'paste':os.path.join(main_script_path,'paste_sums.sh'),'filter':os.path.join(main_script_path,'filter_new_sjs.sh'),'merge':os.path.join(workflow.basedir, 'merge', 'merge.py'),'annotate':os.path.join(workflow.basedir, 'annotate', 'annotate_sjs.py')}

if 'existing_sj_db' not in config:
	config['existing_sj_db']=""
if 'existing_sums' not in config:
	config['existing_sums']=""

wildcard_constraints:
	group_num="\d\d",
	type="(all)|(unique)"

rule all:
	input:
		expand("{file}", file=FILES)

#tar and gzip all the logs, but maintain the directory structure
rule tar_logs:
	input:
		config['input']
	output:
		os.path.join(config['staging'], config['study'] + '.all.logs.tar.gz')
	params:
		log_path=config['input'] + '/' + config['study'] + '/??/*/*.log'
	shell:
		"tar -cvzf {output} {params.log_path}"

#exon sum pasting related rules
rule find_sums:
	input: 
		config['input'], config['sample_ids_file']
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
		config['staging'] + '/{type}.exon_bw_count.{group_num}.decompressed'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging'],
		script_path=SCRIPTS['decompress'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {params.staging}/{params.type}.exon_bw_count.{params.group_num}.manifest {output}"

#do a rule instantiation per low-order name grouping to do hierarchical pastes
rule paste_sums_per_group:
	input:
		config['staging'] + '/{type}.exon_bw_count.{group_num}.decompressed'
	output:
		config['staging'] + '/{type}.exon_bw_count.{group_num}.pasted'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging'],
		script_path=SCRIPTS['paste'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {params.staging}/{params.type}.exon_bw_count.{params.group_num}.manifest {output}"

def get_pasted_sum_files(wildcards):
	return [config['staging']+"/%s.exon_bw_count.%s.pasted" % (wildcards.type, f.split('/').pop()) for f in glob.glob(config['input']+'/*/??')]	

rule collect_pasted_sums:
	input:
		get_pasted_sum_files
	output:
		config['staging'] + '/{type}.exon_bw_count.groups.pasted.files.list'
	params:
		staging=config['staging'],
		type=lambda wildcards: wildcards.type
	shell:
		"rm {params.staging}/*.{params.type}.*.unc && ls {params.staging}/{params.type}.exon_bw_count.*.pasted > {params.staging}/{params.type}.exon_bw_count.groups.pasted.files.list"

rule paste_sums_final:
	input:
		config['staging'] + '/{type}.exon_bw_count.groups.pasted.files.list'
	output:
		os.path.join(config['staging'], config['study'] + '.{type}.exon_bw_count.pasted.gz')
	params:
		staging=config['staging'],
		script_path=SCRIPTS['paste'],
		existing_sums=config['existing_sums'],
		type=lambda wildcards: wildcards.type
	shell:
		"{params.script_path} {input} {output} dont_get_ids {params.existing_sums} && rm {params.staging}/{params.type}.*.pasted"

#junction related rules
rule find_sjs:
	input: 
		config['input'], config['sample_ids_file']
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
		config['staging'] + '/sj.{group_num}.manifest.filtered'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging'],
		script_path=SCRIPTS['filter']
	shell:
		"{params.script_path} {params.staging}/sj.{params.group_num}.manifest"

rule merge_sjs:
	input:
		config['staging'] + '/sj.{group_num}.manifest.filtered'
	output:
		config['staging'] + '/sj.{group_num}.manifest.sj_sample_files.merged.tsv.gz'
	params:
		staging=config['staging'],
		filtered_manifest=lambda wildcards, input: '.'.join(input[0].split('.')[:-1]),
		script_path=SCRIPTS['merge']
	shell:
		"pypy {params.script_path} --list-file {params.filtered_manifest} --gzip | sort -k1,1 -k2,2n -k3,3n | gzip > {params.filtered_manifest}.sj_sample_files.merged.tsv.gz"

def get_sj_merged_files(wildcards):
	return [config['staging']+"/sj.%s.manifest.sj_sample_files.merged.tsv.gz" % f.split('/').pop() for f in glob.glob(config['input']+'/*/??')]	

rule collect_merged_sjs:
	input:
		get_sj_merged_files
	output:
		config['staging'] + '/sj_groups.merged.files.list'
	params:
		staging=config['staging']
	shell:
		"ls {params.staging}/*.sj_sample_files.merged.tsv.gz > {params.staging}/sj_groups.merged.files.list"

rule merge_all_sjs:
	input: 
		config['staging'] + '/sj_groups.merged.files.list'
	output:
		config['staging'] + '/all.sjs.merged.tsv.gz'
	params:
		script_path=SCRIPTS['merge'],
		existing_sj_db=config['existing_sj_db']
	shell:
		"pypy {params.script_path} --list-file {input} --gzip --append-samples --existing-sj-db \"{params.existing_sj_db}\" | gzip > {output}"

rule annotate_all_sjs:
	input:
		config['staging'] + '/all.sjs.merged.tsv.gz'
	output:
		os.path.join(config['staging'], config['study'] + '.all.sjs.merged.annotated.tsv.gz')
	params:
		annot_sjs=config['annotated_sjs'],
		script_path=SCRIPTS['annotate']
	shell:
		"zcat {input} | pypy {params.script_path} --compiled-annotations {params.annot_sjs} | gzip > {output}"
