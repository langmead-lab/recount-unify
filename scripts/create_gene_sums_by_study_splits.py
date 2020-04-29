import sys
import gzip
from datetime import datetime

date_split = str(datetime.today())

compilation = 'sra'
feat_type = 'gene_sums'

idsF = sys.argv[1]
sumsF = sys.argv[2]
#SIRV
annotation = sys.argv[3]
#outdir='supers/gene_sums_per_study/' + annotation
outdir=annotation

run2study = {}
with open(idsF,"r") as fin:
	run2study = {line.split('\t')[1]:line.split('\t')[0] for line in fin}

#posmap = {}
study_cuts = {study:[] for study in run2study.values()}
fin = gzip.open(sumsF,"r")
#get header
fields = fin.readline().rstrip().split('\t')
#posmap = {run:i for (i,run) in enumerate(fields)}
for (i,run) in enumerate(fields):
	if run in run2study:
		study_cuts[run2study[run]].append(i+1)
fin.close()

for study in study_cuts.keys():
	if len(study_cuts[study]) > 0:
		#only cut out the sample counts, we can always add in the gid,chr,bplengthstart,end separately, saves space
		sys.stdout.write("cat <(echo '##annotation=%s') <(echo '##date.generated=%s') <(pigz --stdout -p2 -d %s | cut -f 1,%s) | pigz --fast -p3 > %s/%s.%s.%s.%s.gz\n" % (annotation, date_split, sumsF, (','.join([str(run_pos) for run_pos in study_cuts[study]])), outdir, compilation, feat_type, study, annotation))
	else:
		sys.stderr.write("%s\tno_runs\n" % study)
