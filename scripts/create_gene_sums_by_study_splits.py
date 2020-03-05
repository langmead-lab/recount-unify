import sys
import gzip


idsF = sys.argv[1]
sumsF = sys.argv[2]
#SIRV
annotation = sys.argv[3]
outdir='supers/gene_sums_per_study/' + annotation

run2study = {}
with open(idsF,"r") as fin:
	run2study = {line.split('\t')[1]:line.split('\t')[0] for line in fin}

posmap = {}
study_cuts = {study:[] for study in run2study.values()}
fin = gzip.open(sumsF,"r")
#get header
fields = fin.readline().rstrip().split('\t')
posmap = {run:i for (i,run) in enumerate(fields)}
for (i,run) in enumerate(fields):
	if run in run2study:
		study_cuts[run2study[run]].append(i)
fin.close()

for study in study_cuts.keys():
	if len(study_cuts[study]) > 0:
		sys.stdout.write("zcat %s | cut -f 1-5,%s | pigz --fast -p1 > %s/%s.%s.gz\n" % (sumsF, (','.join([str(run_pos) for run_pos in study_cuts[study]])), outdir, annotation, study))
