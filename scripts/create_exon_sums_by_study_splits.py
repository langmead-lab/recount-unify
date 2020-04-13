import sys
import gzip

#exon version, no individual annotation splits (all kept together for efficiency)
#also, header is separate from sums file

idsF = sys.argv[1]
sumsF = sys.argv[2]
headerF = sys.argv[3]
outdir='supers/exon_sums_per_study/'

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
        #only cut out the sample counts, we can always add in the gid,chr,bplengthstart,end separately, saves space
        cut_fields = (','.join([str(run_pos) for run_pos in study_cuts[study]]))
        sys.stdout.write("cat <(cut -f %s %s) <(pigz --stdout -p 2 -d %s | cut -f %s) | pigz --fast -p1 > %s/%s.gz\n" % (cut_fields, headerF, sumsF, cut_fields, outdir, study))
    else:
        sys.stderr.write("%s\tno_runs\n" % study)
