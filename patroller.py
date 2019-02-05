#!/bin/env python3.6
import sys
import os
import shutil
from time import sleep
import datetime
import glob
from pathlib import Path
import re
import subprocess
import argparse
import urllib.request as urlr
from collections import Counter
import logging

SNAKEMAKE_PATH = 'snakemake'

#in sec
SLEEP=5
LOGS_DIR="logs"
PUMP_DIR="srav1"
DEST_DIR="staging"
OUTPUT_DIR="unified"

#example: https://recount-meta.s3.amazonaws.com/srav1/srav1.txt'
S3_MANIFEST_URL_PREFIX='http://recount-meta.s3.amazonaws.com'

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('patroller')
def log(msg, logger=logger.info):
    logger(msg)

def retrieve_project_manifest(args):
    url = '%s/%s/%s.txt' % (S3_MANIFEST_URL_PREFIX, args.project, args.project)
    log("retrieving study2run manifest from S3: %s" % url)
    with urlr.urlopen(str(url)) as fin:
        study2run_count = Counter(x.rstrip().split(' ')[0] for x in str(fin.read()).split('\\n')[:-1])
        return study2run_count

def run_command(cmd_args, cmd_name, datetime_stamp):
    #dont use shlex.quote, screws up the commandline
    cmd_args = ' '.join([' '.join(cmd_args), '>', '%s/%s.%s' % (LOGS_DIR, cmd_name, datetime_stamp), '2>&1'])
    log("running %s" % cmd_args)
    try:
        cp = subprocess.run(args=cmd_args, shell=True, check=True, universal_newlines=True) 
    except subprocess.CalledProcessError as cpe:
        sys.stderr.write("error in run_command for command: %s\n" % cmd_args)
        raise cpe

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


def process_study(args, study_loworder, study, study_map, sample_ids_file):
    staging_dir = os.path.join(args.staging_dir, study_loworder, study)
    if not os.path.exists(staging_dir):
        Path(staging_dir).mkdir(parents=True, exist_ok=True)
    ds = datetime.datetime.now().timestamp()
    fout = open(os.path.join(staging_dir, "attempts.moved.%s" % (ds)), "w")
    for job in study_map.keys():
        (attempt_num, fpath, run, fname) = study_map[job]
        attempt_name = fpath.split(os.path.sep)[-1]
        run_loworder = run[-2:]
        fout.write("%s\n" % (fpath)) 
        #TODO: this is temporary (making soft links), need to hardlink, then delete in a separate process
        d = os.path.join(staging_dir, run_loworder, attempt_name)
        if not os.path.exists(d):
            Path(d).mkdir(parents=True, exist_ok=True)
        try:
            #os.symlink(fpath, os.path.join(d, f))
            #[os.link(file_, os.path.join(d, file_.split(os.path.sep)[-1])) for file_ in glob.glob('%s/*' % fpath)]
            [os.symlink(file_, os.path.join(d, file_.split(os.path.sep)[-1])) for file_ in glob.glob('%s/*' % fpath)]
        except FileExistsError as fee:
            pass
    fout.close()

    output_dir = os.path.join(args.output_dir, study_loworder, study)
    if not os.path.exists(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    #now fire up snakemake on DEST_DIR/study
    snakemake_cmd_args = generate_snakemake_cmd_args(args, staging_dir, output_dir, study, sample_ids_file)
    run_command(snakemake_cmd_args, study, ds)
    return True


def create_parser():
    parser = argparse.ArgumentParser(description='recount-unify script to "patrol" the directories containing processed files from recount-pump and move them to a pre-staging directorty, and then run the jx merge/annotate and exon sum paste steps (recount-unify) on them.')
    parser.add_argument('--project', metavar='project-name', type=str, required=True, help='name of project (group of one or more studies run together); this is used to retrieve the manifest of study2run mappings to determine which studies are finished and ready to be unified')
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
    parser.add_argument('--debug', action='store_const', const=True, default=False, help='run one time through loop processing up to 5 studies')

    return parser
    
fpath_patt = re.compile(r'^((.+)_attempt(\d+)).done$')
def find_done_runs(args, loworders, studies_done, seen, attempts_tracker, finished_studies):
    files = glob.glob("%s/**/*.done" % (args.incoming_dir), recursive=True)
    count = 0
    for f in files:
        #args.pump_dir/??/study/proj#_input#_attempt#.done
        #old example: srav1/39/SRP008339/proj1_input5472_attempt0.done
        #new example (loworder from run also): geuv_sc/42/ERP001942/92/proj1_input193_attempt0.done
        m = fpath_patt.search(f)
        assert(m)
        #e.g. geuv_sc/42/ERP001942/92/proj1_input193_attempt0
        fdir = m.group(1)
        #e.g. geuv_sc/42/ERP001942/92/proj1_input193
        fkey = m.group(2)
        #e.g. 0
        attempt_num = m.group(3)

        #only get the earliest finished one for stability's sake
        #if fkey in attempts_tracker and attempts_tracker[fkey] < attempt_num:
        #attempts_tracker[fkey] = attempt_num 

        fields = fdir.split(os.path.sep)
        #e.g. proj1_input193_attempt0.done
        fname = fields[-1]
        #e.g. ERP001942
        study = fields[-3]
        if study in finished_studies or (study in seen and fkey in seen[study] and seen[study][fkey][0] < attempt_num):
            continue 
        #e.g. 42
        loworder = fields[-4]
        loworders[study] = loworder

        count += 1

        if study not in seen:
            seen[study]={}
        if fkey not in seen[study] or seen[study][fkey][0] > attempt_num:
            seen[study][fkey] = [attempt_num, fdir, fname, fname]

    runs_done_count = 0
    for study in seen.keys():
        for fkey in seen[study].keys():
            (attempt_num, fdir, run, fname) = seen[study][fkey]
            #now get the run accession
            #need this at this step to be able to further split the
            #input for this study into subgroups based on the low order of the run accession
            for run_manifest in glob.glob("%s/*.manifest" % (fdir)):
                #e.g. SRR306518
                run = run_manifest.split(os.path.sep)[-1].split('_')[0]
                seen[study][fkey][2] = run
                studies_done[study] += 1
                runs_done_count += 1
    return runs_done_count


def load_finished_studies(args):
    #e.g. unified/13/SRP001313/SRP001313.all.sjs.merged.annotated.tsv.gz
    finished_studies = glob.glob("%s/??/*/*.all.sjs.merged.annotated.tsv.gz" % (args.output_dir))
    return set([x.split(os.path.sep)[-2] for x in finished_studies])
  
 
def main():
    parser = create_parser()
    args = parser.parse_args()

    study2run_count = retrieve_project_manifest(args)
    finished_studies = load_finished_studies(args)

    loworders = {}
    seen = {}
    done = {}
    counter = 0
    loop = True
    if args.debug:
        logger.setLevel(logging.DEBUG)
        log("DEBUG mode")
    round_idx = 0
    attempts_tracker = {}
    while(loop):
        runs_by_study_done = Counter()
        log("ROUND %d\tfinding runs which are done in\t%s" % (round_idx, args.incoming_dir))
        runs_done_count = find_done_runs(args, loworders, runs_by_study_done, seen, attempts_tracker, finished_studies)
        log("ROUND %d\truns found which are done this round:\t%d" % (round_idx, runs_done_count))
        #overall, keep track of studies on FS that aren't done vs. those that have been unified already 
        keys = set(seen.keys())
        studies_to_process = keys - finished_studies
        log("ROUND %d\tstudies to process this round:\t%d" % (round_idx, len(studies_to_process)))
        ctr = 0
        while(len(studies_to_process) > 0):
            study = studies_to_process.pop()
            log("ROUND %d\tpopping study\t%s" % (round_idx, study))
            #first check to see if all the runs for this study are done
            #if not, skip for now
            num_runs_done = runs_by_study_done[study]
            num_runs_total = study2run_count[study]
            if num_runs_done != num_runs_total:
                log("ROUND %d\tskipping study\t%s\t%d\tvs\t%d" % (round_idx, study, num_runs_total, num_runs_done))
                continue
            #for debugging
            if args.debug and ctr >= 7:
                break
            ctr += 1
            log("ROUND %d\tprocessing study\t%s" % (round_idx, study))
            #pump processing is done, so prepare the file hierarchy and unify the results
            processed = process_study(args, loworders[study], study, seen[study], args.sample_ID_file) 
            if processed:
                finished_studies.add(study)
                counter += 1
                log("ROUND %d\tprocessed study successfully\t%s" % (round_idx, study))
        #for debugging, only loop once
        loop = False
        #loop = not args.debug
        round_idx += 1
        #log("Sleeping before round %d" % round_idx)
        #sleep(SLEEP)

if __name__ == '__main__':
    main()

