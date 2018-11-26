#!/bin/env python3.6
import sys
import os
import shutil
import datetime
import glob
from pathlib import Path
import re

#in sec
SLEEP=10
PUMP_DIR="./sra"
DEST_DIR="./pre_staging2"

def search_for_files(directory):
    return glob.glob("%s/**/*.done" % (directory), recursive=True)

def check_study_status(study):
    return True

def process_study(study, study_map):
    if not os.path.exists(os.path.join(DEST_DIR, study)):
        Path(os.path.join(DEST_DIR, study)).mkdir()
    fout = open(os.path.join(DEST_DIR, study, "attempts.moved.%s" % (datetime.datetime.now().timestamp())), "w")
    for job in study_map.keys():
        (fpath, loworder, attempt_num) = study_map[job]
        #if not os.path.exists(os.path.join(DEST_DIR, study, loworder)):
        #    Path(os.path.join(DEST_DIR, study, loworder)).mkdir()
        fout.write("%s\n" % (fpath)) 
        shutil.copytree(fpath, os.path.join(DEST_DIR, study, loworder, job+'_'+attempt_num ))
        #shutil.rmtree(fpath)
    fout.close()

    #now fire up snakemake on DEST_DIR/study
        
def main():
    seen = {}
    fpath_patt = re.compile(r'^((.+)_attempt(\d+)).done$')
    loop = True
    while(loop):
        files = search_for_files(PUMP_DIR)
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
            seen[study][fkey] = [fdir, loworder, attempt_num]
        
        for study in seen.keys():
            status = check_study_status(study)
            #not done processing through the pump
            if not status:
                continue
            #pump processing is done, so prepare the file hierarchy and unify the results
            process_study(study, seen[study]) 
            del seen[study]

        loop = False

if __name__ == '__main__':
    main()

