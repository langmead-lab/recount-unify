#!/bin/env python3.6
import sys
import os
import shutil
import datetime
import glob
from pathlib import Path
import re
import subprocess

#stubbing the basic set of compilations here before we code the
#interface to the remote DB
compilations={'sra':[1,'sra_samples.tsv',1],
                'gtex':[2,'gtex_samples.tsv',1],
                'tcga':[3,'gtex_samples.tsv',23]}

SNAKEMAKE_PATH = 'snakemake'

#in sec
SLEEP=10
PUMP_DIR="./sra"
DEST_DIR="./pre_staging"

def run_command(cmd_args, cmd_name, datetime_stamp):
    for param in ['>', '%s.%s' % (cmd_name, datetime_stamp), '2>&1']:
        cmd_args.append(param)
    #cmd_args.append(" > %s.%s 2>&1" % (cmd_name, datetime_stamp))
    print ("running %s" % cmd_args.join(" "))
    try:
        cp = subprocess.run(cmd_args, shell=True, check=True, text=True) 
    except CalledProcesserror cpe:
        sys.stdout.write("error in run_command for command: %s\n" % (cmd_args.join(" ")))
        raise cpe

def search_for_files(directory):
    return glob.glob("%s/**/*.done" % (directory), recursive=True)

#stubbed currently, needs to interface with remote DB in AWS
#needs to know:
#0) is study fully done processing?
#1) study2compilation mapping
#2) compilation ID prefix code
#3) compilation last used counts (per loworder prefix)
#
#from this set of data it can either:
#1) use existing study sample IDs if present
#2) generate study sample IDs if not present already
#this assumes we only ever process one study at a time at the recount-unify stage
def check_study_status(args, study):
    sample_ids_file = args.sample_ID_file
    if sample_ids_file is None:
        sample_ids_file = os.path.join(dest_dir, "%s.samples.tsv" % (study))
        if not os.path.exists(sample_ids_file):
            #generate sample ids file for study
            pass
    return (True, sample_ids_file)

def generate_snakemake_cmd_args(args, output_dir, sample_ids_file):
    cmd_args = [SNAKEMAKE_PATH,'-j%d' % args.num_sm_proc, args.snakefile, '-r', '--config']
    cmd_args.append('input=%s' % args.pre_staging_dir)
    cmd_args.append('staging=%s' % output_dir)
    cmd_args.append('sample_ids_file=%s' % sample_ids_file)
    cmd_args.append('annotated_sjs=%s' % args.annotated_sj_file)

    #optional config params
    if args.existing_sj_file is not None:
        cmd_args.append('existing_sj_db=%s' % args.existing_sj_file)
    if args.existing_exon_sums_file is not None:
        cmd_args.append('existing_sums=%s' % args.existing_exon_sums_file)
    return cmd_args


def process_study(args, study, study_map, sample_ids_file):
    dest_dir = args.pre_staging_dir
    if not os.path.exists(os.path.join(dest_dir, study)):
        Path(os.path.join(dest_dir, study)).mkdir()
    ds = datetime.datetime.now().timestamp()
    fout = open(os.path.join(dest_dir, study, "attempts.moved.%s" % (ds)), "w")
    for job in study_map.keys():
        (attempt_num, fpath, loworder) = study_map[job]
        fout.write("%s\n" % (fpath)) 
        shutil.copytree(fpath, os.path.join(dest_dir, study, loworder, job+'_'+attempt_num ))
        if args.remove_incoming:
            shutil.rmtree(fpath)
    fout.close()

    output_dir = os.path.join(args.output_dir, study)
    if not os.path.exists(output_dir):
        Path(output_dir).mkdir()
    #now fire up snakemake on DEST_DIR/study
    snakemake_cmd_args = generate_snakemake_cmd_args(args, output_dir, sample_ids_file)
    run_command(snakemake_cmd_args, 'snakemake', ds)

        
def main():
    parser = argparse.ArgumentParser(description='recount-unify script to "patrol" the directories containing processed files from recount-pump and move them to a pre-staging directorty, and then run the jx merge/annotate and exon sum paste steps (recount-unify) on them.')
    parser.add_argument('--incoming-dir', metavar='/path/to/dir_containing_pump_processed_files', type=str, default=PUMP_DIR, help='the path where recount-pump dumps it\'s procssed files.')
    parser.add_argument('--pre-staging-dir', metavar='/path/to/dir_to_store_finished_pump_processed_files', type=str, default=DEST_DIR, help='the path where recount-pump runs which have the *.done file present, are moved.')
    parser.add_argument('--output-dir', metavar='/path/to/dir_to_store_recount-unify_output', type=str, default=DEST_DIR, help='the path where in process and fully merged sj and exon sum files are stored.')
    parser.add_argument('--remove-incoming', action='store_const', const=True, default=False, help='removes the source files from the recount-pump incoming directories after they\'ve been copied to the pre-staging area')
    parser.add_argument('--annotated-sj-file', metavar='/path/to/file_with_annotated_junctions', type=str, required=True, help='this file contains a list of junctions from one or more annotations.')
    parser.add_argument('--existing-sj-file', metavar='/path/to/file_with_existing_junctions', type=str, default=None, help='if merging these new splice junctions with an already merged set of splice junctions')
    parser.add_argument('--existing-exon-sums-file', metavar='/path/to/file_with_existing_exon_sums', type=str, default=None, help='if merging these new exon sums with an already pasted set of exon sample sums OR just the original BED file/coordinates.')
    parser.add_argument('--sample-ID-file', metavar='/path/to/file_with_sample_IDs_mapped_to_study_and_run_accessions', type=str, default=None, help='this file contains the pre-generated mapping between recount-unify assigned sample IDs and study/run accessions/UUIDs. If not provided and not already present in the default location (<pre-staging>/<study-accession>.samples.tsv) it will be generated.')
    parser.add_argument('--snakefile', metavar='/path/to/file_with_snakemake_rules', type=str, required=True, help='this file contains the list of rules making up recount-unify (sj merging/annotating & exon sum pasting).')
    parser.add_argument('--num-sm-proc', metavar='1', type=int, default=1, help='# of Snakemake parallel processes to run, defaults to 1')

    seen = {}
    fpath_patt = re.compile(r'^((.+)_attempt(\d+)).done$')
    loop = True
    while(loop):
        files = search_for_files(args.pump_dir)
        for f in files:
            #proj1_input44_attempt7.done
            m = fpath_patt.search(f)
            assert(m)
            fdir = m.group(1)
            fkey = m.group(2)
            attempt_num = m.group(3)

            fields = f.split('/')
            fname = fields[-1]
            loworder = fields[-2]
            study = fields[-3]

            if study not in seen:
                seen[study]={}
            if fkey not in seen[study] or seen[study][fkey][0] > attempt_num:
                seen[study][fkey] = [attempt_num, fdir, loworder]
        
        for study in seen.keys():
            #find out:
            #1) study started processing?
            #2) have sample IDs been generated for the study, and if so where are they?
            #if no to 2), this will generate the sample IDs for the whole study
            (status, sample_ids_file) = check_study_status(args, study)
            #not done processing through the pump
            if not status:
                continue
            #pump processing is done, so prepare the file hierarchy and unify the results
            process_study(args, study, seen[study], sample_ids_file) 
            del seen[study]

        loop = False

if __name__ == '__main__':
    main()

