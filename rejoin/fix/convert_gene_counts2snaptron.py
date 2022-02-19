#!/usr/bin/env python3
import sys

#does 2 things:
#1) produced a snaptron-compatbile list of a) sample IDs (rail ids) and b) counts for all genes which have any coverage

#NO HEADER
study_samplesF=sys.argv[1]
id_mapF=sys.argv[2]
#report_outputF=sys.argv[3]

sample_positions={}
with open(study_samplesF,"r") as fin:
    sample_positions={sample.rstrip():i for i,sample in enumerate(fin)}

#report_out = open(report_outputF,"w")

id_map = {}
with open(id_mapF,"r") as fin:
    for line in fin:
        line = line.rstrip()
        (rail_id,external_id,study_id)=line.split('\t')
        #may be samples in the samples.tsv which aren't actually in recount3 counts
        if external_id in sample_positions:
            id_map[sample_positions[external_id]] = rail_id

for line in sys.stdin:
    line_original = line
    line = line.rstrip()
    fields = line.split('\t')
    snaptron_counts = [id_map[i]+':'+count for i,count in enumerate(fields[1:]) if int(count) > 0]
    if len(snaptron_counts) > 0:
        sys.stdout.write(fields[0]+'\t'+','+','.join(snaptron_counts)+'\n')
    else:
        #output anyway as we're going to do a paste across studies and study groups
        sys.stdout.write(fields[0]+'\t\n')
