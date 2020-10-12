import sys
import gzip
from datetime import datetime

date_split = str(datetime.today())

#annotation='M023'

#exon version, no individual annotation splits (all kept together for efficiency for now)
#also, header is separate from sums file

idsF = sys.argv[1]
sumsF = sys.argv[2]
headerF = sys.argv[3]
#outdir='supers/exon_sums_per_study/'
outdir = sys.argv[4]

run2study = {}
with open(idsF,"r") as fin:
    run2study = {line.split('\t')[1]:line.split('\t')[0] for line in fin}

study_cuts = {study:[] for study in run2study.values()}

with open(headerF,"r") as fin:
    #get header
    fields = fin.readline().rstrip().split('\t')
    for (i,run) in enumerate(fields):
        if run in run2study:
            study_cuts[run2study[run]].append(i+1)

for study in study_cuts.keys():
    if len(study_cuts[study]) > 0:
        #find contiguous ranges of 
        #sample positions to condense
        #this keeps the command line smaller
        prev_run_pos = -100
        start_pos = -100
        positions = []
        for run_pos in study_cuts[study]:
            if run_pos != prev_run_pos+1:
                pos = prev_run_pos
                if start_pos != pos:
                    pos = "%d-%d" % (start_pos, prev_run_pos)
                if prev_run_pos != -100:
                    positions.append(str(pos))
                start_pos = run_pos 
            prev_run_pos = run_pos
        pos = prev_run_pos
        if start_pos != pos:
            pos = "%d-%d" % (start_pos, prev_run_pos)
        if prev_run_pos != -100:
            positions.append(str(pos))
        cut_fields = ','.join(positions) 
        #only cut out the sample counts, we can always add in the gid,chr,bplengthstart,end separately, saves space
        #cut_fields = (','.join([str(run_pos) for run_pos in study_cuts[study]]))
        #since we're not doing the annotation split here don't include the extra headers yet
        sys.stdout.write("cat <(cut -f 1,%s %s) <(pigz --stdout -p 2 -d %s | cut -f 1,%s) | pigz --fast -p2 > %s/%s.gz\n" % (cut_fields, headerF, sumsF, cut_fields, outdir, study))
    else:
        sys.stderr.write("%s\tno_runs\n" % study)
