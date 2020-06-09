#!/usr/bin/env python
"""
annotate_jxs.py
Based on Abhi Nellore's script process_introns.py.
Now just annotate intron splice sites (no sequence)
and include the length of the intron and compilation
specific snaptron_id (assigned here).

Reads introns from stdin (all_SRA_introns.tsv.gz), 
annotation files from command-line parameters and annotation from
all_SRA_introns.tsv.gz should have AT LEAST the following tab-separated fields
on each line:
0. chromosome
1. start position
2. end position
3. strand
4. 5' motif (left or donor)
5. 3' motif (right or acceptor)
(can include other fields after this, but are simply passed through
and appended to the end of the output line)

Tab-separated output written to stdout (unchanged unless noted):
0. id
1. chromosome
2. start position
3. end position
4. length of intron
5. strand
6. 1 if junction is annotated else 0
7. left splice seq
8. right splice seq
9. 1 if left splice site is annotated else 0
10. 1 if right splice site is annotated else 0

everything else that was after the 3' motif column

additionally, sample coverage related fields:
count of samples
sum
average
median

compilation_id (data_source_id)
"""
import sys
import os
import subprocess
import struct
from operator import itemgetter
from bisect import bisect_right
from collections import defaultdict
import string
import gzip

CHRM_COL=0
START_COL=1
END_COL=2
STRAND_COL=3
MOTIF_COL=4
SAMP_COV_COL=5

def load_preformatted_annotated_junctions(f):
    annotated_junctions = {}
    five_p = {}
    three_p = {}
    with gzip.open(f,"r") as fin:
        for line in fin:
            (chr_,start,end,strand,sources) =(None,None,None,None,None)
            fields_ = line.rstrip().split("\t")
            if(len(fields_) > STRAND_COL+1): 
                (chr_,start,end,strand,sources) = line.rstrip().split("\t")
            else:
                (chr_,start,end,strand) = line.rstrip().split("\t")
                sources = "all"
            hkey = (chr_,start,end)
            sources_set = set(sources.split(","))
            if hkey not in annotated_junctions:
                annotated_junctions[hkey]=set()
            annotated_junctions[hkey].update(sources_set)
            hkey = (chr_,start)
            if hkey not in five_p:
                five_p[hkey]=set()
            five_p[hkey].update(sources_set)
            hkey = (chr_,end)
            if hkey not in three_p:
                three_p[hkey]=set()
            three_p[hkey].update(sources_set)
    return (annotated_junctions, five_p, three_p)

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
             
if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--compiled-annotations', type=str, required=True,
            help='TSV of pre-compiled junctions which occur in >=1 annotations'
        )
    parser.add_argument('--compilation-id', type=int, required=True,
            help='source compilation ID for Snaptron output'
        )
    args = parser.parse_args()

    (annotated_junctions, five_p, three_p) = load_preformatted_annotated_junctions(args.compiled_annotations) 

    snaptron_id = 0
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        #allows us to skip the intermediate rearrangement script
        #assuming we're getting this from the output of perbase in motif mode
        tokens[MOTIF_COL] = tokens.pop().upper()
        junction = tuple(tokens[:STRAND_COL])
        annotated = set()
        if junction in annotated_junctions:
            annotated = set(["1"])
        chrm = junction[CHRM_COL]
        start = junction[START_COL]
        end = junction[END_COL]
        length = (int(end) + 1) - int(start)
        strand = tokens[STRAND_COL]
        #sys.stderr.write(line)
        (left_motif, right_motif) = tokens[MOTIF_COL].split('-')
        additional_fields = tokens[SAMP_COV_COL:]

        if (chrm, start) in five_p:
            five_atype = five_p[(chrm, start)]
        if (chrm, end) in three_p:
            three_atype = three_p[(chrm, end)]
            
        covs = [int(s.split(':')[1]) for s in tokens[SAMP_COV_COL].split(',')[1:]]
        sum_ = sum(covs)
        (count, avg, median) = sample_summary_stats(sum_, covs)

        additional_fields.extend([str(count), str(sum_), str(avg), str(median), str(args.compilation_id)])

        print '\t'.join([str(snaptron_id), '\t'.join(junction), str(length), strand,
            ",".join(sorted(annotated)) if len(annotated) > 0 else '0', left_motif, right_motif,
            "%s" % (",".join(sorted(five_atype))) if (chrm, start) in five_p else '0',
            "%s" % (",".join(sorted(three_atype))) if (chrm, end) in three_p else '0',
            '\t'.join(additional_fields)])

        snaptron_id+=1
