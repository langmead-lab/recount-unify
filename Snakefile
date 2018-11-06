#start
FILES=[config['staging'] + '/all.jxs.merged.annotated.tsv.gz']

rule all:
	input:
		expand("{file}", file=FILES)

rule find_jxs:
	input: 
		config['input'], config['sample_ids_file']
	output:
		config['staging'] + '/groups.manifest'
	params:
		staging=config['staging']
	shell:
		"./find_new_jxs.sh {input[0]} {input[1]} {params.staging}"

rule filter_jxs:
	input:
		config['staging'] + '/groups.manifest'
	output:
		config['staging'] + '/{group_num}.manifest.filtered'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging']
	shell:
		"./filter_new_jxs.sh {params.staging}/{params.group_num}.manifest"

rule merge_jxs:
	input:
		config['staging'] + '/{group_num}.manifest.filtered'
	output:
		config['staging'] + '/{group_num}.manifest.jx_sample_files.merged.tsv.gz'
	params:
		staging=config['staging'],
		filtered_manifest=lambda wildcards, input: '.'.join(input[0].split('.')[:-1])
	shell:
		"python2 merge.py --list-file {params.filtered_manifest} --gzip | gzip > {params.filtered_manifest}.jx_sample_files.merged.tsv.gz"

import glob
def get_jx_merged_files(wildcards):
	return [config['staging']+"/%s.manifest.jx_sample_files.merged.tsv.gz" % f.split('/').pop() for f in glob.glob(config['input']+'/*/??')]	

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
