#start
#FILES=[config['staging'] + '/all.sjs.merged.annotated.tsv.gz', config['staging'] + '/sums.all.pasted.gz']
FILES=[config['staging'] + '/all.sjs.merged.annotated.tsv.gz'] #, config['staging'] + '/sums.all.pasted.gz']
main_script_path=os.path.join(workflow.basedir,'scripts')
SCRIPTS={'find':os.path.join(main_script_path,'find_new_files.sh'),'paste':os.path.join(main_script_path,'paste_sums.sh'),'filter':os.path.join(main_script_path,'filter_new_sjs.sh'),'merge':os.path.join(workflow.basedir, 'merge', 'merge.py'),'annotate':os.path.join(workflow.basedir, 'annotate', 'annotate_sjs.py')}

if 'existing_sj_db' not in config:
	config['existing_sj_db']=""
if 'existing_sums' not in config:
	config['existing_sums']=""

rule all:
	input:
		expand("{file}", file=FILES)

rule find_sums:
	input: 
		config['input'], config['sample_ids_file']
	output:
		config['staging'] + '/sums.groups.manifest'
	params:
		staging=config['staging'],
		script_path=SCRIPTS['find']
	shell:
		"{params.script_path} {input[0]} {input[1]} {params.staging} sums"

#do a rule instantiation per low-order name grouping to do hierarchical pastes
rule paste_sums_per_group:
	input:
		config['staging'] + '/sums.groups.manifest'
	output:
		config['staging'] + '/sums.{group_num}.pasted'
	params:
		group_num=lambda wildcards: wildcards.group_num,
		staging=config['staging'],
		script_path=SCRIPTS['paste']
	shell:
		"{params.script_path} {params.staging}/sums.{params.group_num}.manifest {output}"

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
		config['staging'] + '/sums.all.pasted.gz'
	params:
		script_path=SCRIPTS['paste'],
		existing_sums=config['existing_sums']
	shell:
		"{params.script_path} {input} {output} dont_get_ids {params.existing_sums}"
			

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
		config['staging'] + '/all.sjs.merged.annotated.tsv.gz'
	params:
		annot_sjs=config['annotated_sjs'],
		script_path=SCRIPTS['annotate']
	shell:
		"zcat {input} | pypy {params.script_path} --compiled-annotations {params.annot_sjs} | gzip > {output}"
