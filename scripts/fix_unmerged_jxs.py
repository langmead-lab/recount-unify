#!/usr/bin/env python
import sys

def sample_summary_stats(sum_, covs):
    covs = sorted(covs)
    count = len(covs)
    median = int(count/2)
    if count % 2 == 0:
        median = round((covs[median-1] + covs[median])/2.0, 3)
    else:
        median = covs[median]
    avg = round(sum_/float(count), 3)
    return (count, avg, median)

SAMP_COV_COL = 11

count_col = SAMP_COV_COL + 1
sum_col = SAMP_COV_COL + 2

prev_key = None
prev_line = None
prev_samps = ""
prev_snid = None
running = 0
cid = None
#chr1    11527   184862  173336  -       0       CT      AC      0       0       ,12:2	11      67      6.091   5       2
for line in sys.stdin:
    fields = line.rstrip().split('\t')
    key = '\t'.join(fields[1:SAMP_COV_COL])
    #fields[count_col] = int(fields[count_col])
    #fields[sum_col] = int(fields[sum_col])
    if key == prev_key:
        #pcount += fields[count_col]
        #psum += fields[sum_col]
        prev_samps += fields[SAMP_COV_COL]
        running += 1
        continue
    #didn't match, print out
    if prev_key is not None:
        if running > 0:
            covs = [int(s.split(':')[1]) for s in prev_samps.split(',')[1:]]
            sum_ = sum(covs)
            (count, avg, median) = sample_summary_stats(sum_, covs)
            cid = fields[-1]
            sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (prev_snid, prev_key, prev_samps, str(count), str(sum_), str(avg), str(median), cid))
        else:
            sys.stdout.write(prev_line)

    #now set the new key
    running = 0
    prev_samps = ""
    prev_key = key
    prev_line = line
    prev_snid = fields[0]
    #pcount = fields[sum_col]
    #psum = fields[sum_col]

if prev_key is not None:
    if running > 0:
        covs = [int(s.split(':')[1]) for s in prev_samps.split(',')[1:]]
        sum_ = sum(covs)
        (count, avg, median) = sample_summary_stats(sum_, covs)
        cid = fields[-1]
        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (prev_snid, prev_key, prev_samps, str(count), str(sum_), str(avg), str(median), cid))
    else:
        sys.stdout.write(prev_line)
