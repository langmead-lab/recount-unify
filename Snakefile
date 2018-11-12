#start
FILES=[config['staging'] + '/all.jxs.merged.annotated.tsv.gz', config['staging'] + '/sums.all.pasted']

rule all:
	input:
		expand("{file}", file=FILES)


rule find_sums:
	input: 
		config['input'], config['sample_ids_file']
	output:
		config['staging'] + '/sums.groups.manifest'
	params:
		staging=config['staging']
	shell:
		"./find_new_files.sh {input[0]} {input[1]} {params.staging} sums"

#do a rule instantiation per low-order name grouping to do hierarchical pastes
rule paste_sums_per_group:
	input:
		config['staging'] + '/sums.groups.manifest'
	output:
		config['staging'] + '/sums.{group_num}.pasted'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging']
	shell:
		"./paste_sums.sh {params.staging}/sums.{params.group_num}.manifest {output}"

import glob
def get_pasted_sum_files(wildcards):
	return [config['staging']+"/sums.%s.pasted" % f.split('/').pop() for f in glob.glob(config['input']+'/*/??')]	

rule collect_pasted_sums:
	input:
		get_pasted_sum_files
	output:
		config['staging'] + '/sums_groups.pasted.files.list'
	params:
		staging=config['staging']
	shell:
		"ls {params.staging}/sums.*.pasted > {params.staging}/sums_groups.pasted.files.list"

rule paste_sums_final:
	input:
		config['staging'] + '/sums_groups.pasted.files.list'
	output:
		config['staging'] + '/sums.all.pasted'
	shell:
		"./paste_sums.sh {input} {output}"
			

#junction related rules
rule find_jxs:
	input: 
		config['input'], config['sample_ids_file']
	output:
		config['staging'] + '/jx.groups.manifest'
	params:
		staging=config['staging'],
		wildc='"*.gz"'
	shell:
		"./find_new_files.sh {input[0]} {input[1]} {params.staging} jx {params.wildc}"

rule filter_jxs:
	input:
		config['staging'] + '/jx.groups.manifest'
	output:
		config['staging'] + '/jx.{group_num}.manifest.filtered'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging']
	shell:
		"./filter_new_jxs.sh {params.staging}/jx.{params.group_num}.manifest"

rule merge_jxs:
	input:
		config['staging'] + '/jx.{group_num}.manifest.filtered'
	output:
		config['staging'] + '/jx.{group_num}.manifest.jx_sample_files.merged.tsv.gz'
	params:
		staging=config['staging'],
		filtered_manifest=lambda wildcards, input: '.'.join(input[0].split('.')[:-1])
	shell:
		"python2 merge.py --list-file {params.filtered_manifest} --gzip | gzip > {params.filtered_manifest}.jx_sample_files.merged.tsv.gz"

def get_jx_merged_files(wildcards):
	return [config['staging']+"/jx.%s.manifest.jx_sample_files.merged.tsv.gz" % f.split('/').pop() for f in glob.glob(config['input']+'/*/??')]	

rule collect_merged_jxs:
	input:
		get_jx_merged_files
	output:
		config['staging'] + '/jx_groups.merged.files.list'
	params:
		staging=config['staging']
	shell:
		"ls {params.staging}/*.jx_sample_files.merged.tsv.gz > {params.staging}/jx_groups.merged.files.list"

rule merge_all_jxs:
	input: 
		config['staging'] + '/jx_groups.merged.files.list'
	output:
		config['staging'] + '/all.jxs.merged.tsv.gz'
	shell:
		"python2 ./merge.py --list-file {input} --gzip --append-samples | gzip > {output}"

rule annotate_all_jxs:
	input:
		config['staging'] + '/all.jxs.merged.tsv.gz'
	output:
		config['staging'] + '/all.jxs.merged.annotated.tsv.gz'
	params:
		annot_jxs=config['annotated_jxs']
	shell:
		"zcat {input} | python2 ./annotate_jxs.py --compiled-annotations {params.annot_jxs} | gzip > {output}"
