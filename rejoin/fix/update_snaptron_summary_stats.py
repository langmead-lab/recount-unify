#!/usr/bin/env python3
import sys
import numpy as np

gene2snaptron_samplesF = sys.argv[1]
#1 each for srav3h (0), gtexv2 (2), tcgav2 (3)
dataset_id=sys.argv[2]

#updating: 1) samples_count, 2) coverage_sum, 3) coverage_avg, 4) coverage_median
fin = open(gene2snaptron_samplesF,"r")
for line in fin:
    line = line.rstrip()
    (gid, samples) = line.split('\t')
    sfields = samples.split(',')
    
    #based on snaptron.snaputil:extract_sids_and_covs_from_search_iter() 
    found = np.empty((0,2),dtype=np.int32)
    for sample in sfields[1:]:
        (sid, cov)=sample.split(':')
        found = np.append(found, [[int(sid),int(cov)]], 0)
    length = len(found)
    ssum = np.sum(found[0:length,1])
    savg = np.mean(found[0:length,1])
    smedian = np.median(found[0:length,1])
    sys.stdout.write(line + '\t' + str(length) + '\t' + str(ssum) + '\t' + str(savg) + '\t' + str(smedian) + '\t' + dataset_id + '\n')
