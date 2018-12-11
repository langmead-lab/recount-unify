#!/bin/env python3.6
import sys
import os
import shutil
import datetime
import glob
from pathlib import Path
import re
import subprocess
import argparse
import shlex

#stubbing the basic set of compilations here before we code the
#interface to the remote DB
compilations={'sra':[1,'sra_samples.tsv',1],
                'gtex':[2,'gtex_samples.tsv',1],
                'tcga':[3,'gtex_samples.tsv',23]}

SNAKEMAKE_PATH = 'snakemake'

#in sec
SLEEP=10
LOGS_DIR="logs"
PUMP_DIR="srav1"
DEST_DIR="staging"
OUTPUT_DIR="unified"




def run_command(cmd_args, cmd_name, datetime_stamp):
    #dont use shlex.quote, screws up the commandline
    cmd_args = ' '.join([' '.join(cmd_args), '>', '%s/%s.%s' % (LOGS_DIR, cmd_name, datetime_stamp), '2>&1'])
    print ("running %s" % cmd_args)
    try:
        cp = subprocess.run(args=cmd_args, shell=True, check=True, universal_newlines=True) 
    except subprocess.CalledProcessError as cpe:
        sys.stderr.write("error in run_command for command: %s\n" % cmd_args)
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
    #TODO: return false for now
    return (True, sample_ids_file)

def generate_snakemake_cmd_args(args, staging_dir, output_dir, study, sample_ids_file):
    cmd_args = [SNAKEMAKE_PATH,'-j%d' % args.num_sm_proc, '--snakefile', args.snakefile, '-r', '--config']
    cmd_args.append('study=%s' % study)
    cmd_args.append('input=%s' % staging_dir)
    cmd_args.append('staging=%s' % output_dir)
    cmd_args.append('sample_ids_file=%s' % sample_ids_file)
    cmd_args.append('annotated_sjs=%s' % args.annotated_sj_file)

    #optional config params
    if args.existing_sj_file is not None:
        cmd_args.append('existing_sj_db=%s' % args.existing_sj_file)
    if args.existing_exon_sums_file is not None:
        cmd_args.append('existing_sums=%s' % args.existing_exon_sums_file)
    return cmd_args


def process_study(args, loworder, study, study_map, sample_ids_file):
    staging_dir = os.path.join(args.staging_dir, loworder, study)
    if not os.path.exists(staging_dir):
        Path(staging_dir).mkdir(parents=True, exist_ok=True)
    ds = datetime.datetime.now().timestamp()
    fout = open(os.path.join(staging_dir, "attempts.moved.%s" % (ds)), "w")
    for job in study_map.keys():
        (attempt_num, fpath, run, fname) = study_map[job]
        f = fpath.split(os.path.sep)[-1]
        loworder = run[-2:]
        fout.write("%s\n" % (fpath)) 
        #TODO: this is temporary (making soft links), need to hardlink, then delete
        d = os.path.join(staging_dir, loworder)
        if not os.path.exists(d):
            Path(d).mkdir(parents=True, exist_ok=True)
        try:
            os.symlink(fpath, os.path.join(d, f))
        except FileExistsError as fee:
            pass
        #shutil.copytree(fpath, os.path.join(dest_dir, loworder, study, job+'_'+attempt_num ))
        #if args.remove_incoming:
        #    shutil.rmtree(fpath)
    fout.close()

    output_dir = os.path.join(args.output_dir, loworder, study)
    if not os.path.exists(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    #now fire up snakemake on DEST_DIR/study
    snakemake_cmd_args = generate_snakemake_cmd_args(args, staging_dir, output_dir, study, sample_ids_file)
    print (snakemake_cmd_args)
    run_command(snakemake_cmd_args, study, ds)
    return True

        
def main():
    parser = argparse.ArgumentParser(description='recount-unify script to "patrol" the directories containing processed files from recount-pump and move them to a pre-staging directorty, and then run the jx merge/annotate and exon sum paste steps (recount-unify) on them.')
    parser.add_argument('--incoming-dir', metavar='/path/to/dir_containing_pump_processed_files', type=str, default=PUMP_DIR, help='the path where recount-pump dumps it\'s procssed files.')
    parser.add_argument('--staging-dir', metavar='/path/to/dir_to_store_finished_pump_processed_files', type=str, default=DEST_DIR, help='the path where recount-pump runs which have the *.done file present, are moved.')
    parser.add_argument('--output-dir', metavar='/path/to/dir_to_store_recount-unify_output', type=str, default=OUTPUT_DIR, help='the path where in process and fully merged sj and exon sum files are stored.')
    parser.add_argument('--remove-incoming', action='store_const', const=True, default=False, help='removes the source files from the recount-pump incoming directories after they\'ve been copied to the pre-staging area')
    parser.add_argument('--snakefile', metavar='/path/to/file_with_snakemake_rules', type=str, required=True, help='this file contains the list of rules making up recount-unify (sj merging/annotating & exon sum pasting).')
    parser.add_argument('--num-sm-proc', metavar='1', type=int, default=1, help='# of Snakemake parallel processes to run, defaults to 1')
    #snakemake specific parameters
    parser.add_argument('--annotated-sj-file', metavar='/path/to/file_with_annotated_junctions', type=str, required=True, help='this file contains a list of junctions from one or more annotations.')
    parser.add_argument('--existing-sj-file', metavar='/path/to/file_with_existing_junctions', type=str, default=None, help='if merging these new splice junctions with an already merged set of splice junctions')
    parser.add_argument('--existing-exon-sums-file', metavar='/path/to/file_with_existing_exon_sums', type=str, default=None, help='if merging these new exon sums with an already pasted set of exon sample sums OR just the original BED file/coordinates.')
    parser.add_argument('--sample-ID-file', metavar='/path/to/file_with_sample_IDs_mapped_to_study_and_run_accessions', type=str, default=None, help='this file contains the pre-generated mapping between recount-unify assigned sample IDs and study/run accessions/UUIDs. If not provided and not already present in the default location (<pre-staging>/<study-accession>.samples.tsv) it will be generated.')

    args = parser.parse_args()

    seen = {}
    done = {}
    fpath_patt = re.compile(r'^((.+)_attempt(\d+)).done$')
    loop = True
    counter = 0
    loworders = {}
    while(loop):
        files = search_for_files(args.incoming_dir)
        for f in files:
            #args.pump_dir/??/study/proj#_input#_attempt#.done
            #example: srav1/39/SRP008339/proj1_input5472_attempt0.done
            m = fpath_patt.search(f)
            assert(m)
            #e.g. srav1/39/SRP008339/proj1_input5472_attempt0
            fdir = m.group(1)
            #e.g. srav1/39/SRP008339/proj1_input5472
            fkey = m.group(2)
            #e.g. 1
            attempt_num = m.group(3)
            
            #now get the run accession
            #need this at this step to be able to further split the
            #input for this study into subgroups based on the low order of the run accession
            manifest = glob.glob("%s/*.manifest" % (fdir))[0]
            #e.g. SRR306518
            run = manifest.split(os.path.sep)[-1].split('_')[0]

            fields = fdir.split(os.path.sep)
            #e.g. proj1_input5472_attempt0.done 
            fname = fields[-1]
            #e.g. SRP008339
            study = fields[-2]
            #e.g. 39
            loworder = fields[-3]
            loworders[study] = loworder

            if study not in seen:
                seen[study]={}
            if fkey not in seen[study] or seen[study][fkey][0] > attempt_num:
                seen[study][fkey] = [attempt_num, fdir, run, fname]
       
        keys = set(seen.keys())
        to_process = keys - set(done.keys())
        ctr = 0
        while(len(to_process) > 0):
            study = to_process.pop()
            #find out:
            #1) study started processing?
            #2) have sample IDs been generated for the study, and if so where are they?
            #if no to 2), this will generate the sample IDs for the whole study
            (status, sample_ids_file) = check_study_status(args, study)
            #not done processing through the pump
            if not status or ctr >= 5:
            #if not status:
                continue
            ctr += 1
            #pump processing is done, so prepare the file hierarchy and unify the results
            processed = process_study(args, loworders[study], study, seen[study], sample_ids_file) 
            if processed:
                done[study] = counter
                counter += 1

        loop = False

if __name__ == '__main__':
    main()

