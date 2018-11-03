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

def load_preformatted_annotated_junctions(f):
    annotated_junctions = {}
    five_p = {}
    three_p = {}
    with gzip.open(f,"r") as fin:
        for line in fin:
            (chr_,start,end,strand,sources) =(None,None,None,None,None)
            fields_ = line.rstrip().split("\t")
            if(len(fields_) >= 5): 
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
             
if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--compiled-annotations', type=str, required=True,
            help='TSV of pre-compiled junctions which occur in >=1 annotations'
        )
    args = parser.parse_args()

    (annotated_junctions, five_p, three_p) = load_preformatted_annotated_junctions(args.annotations[0]) 

    snaptron_id = 0
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        junction = tuple(tokens[:3])
        #check to see if we want this junction
        annotated = set()
        if junction in annotated_junctions:
            annotated = set(["1"])
        start = int(junction[1])
        length = int(junction[2]) + 1 - start
        strand = tokens[3]
        tokens_length = len(tokens)
        left_motif, right_motif = tokens[4], tokens[5]
        additional_fields = tokens[6:]
        if (junction[0], junction[1]) in five_p:
            five_atype = five_p[(junction[0], junction[1])]
        if (junction[0], junction[2]) in three_p:
            three_atype = three_p[(junction[0], junction[2])]
        print '\t'.join([str(snaptron_id), '\t'.join(junction), str(length), strand,
            ",".join(sorted(annotated)) if len(annotated) > 0 else '0', left_motif, right_motif,
            "%s" % (",".join(sorted(five_atype))) if (junction[0], junction[1]) in five_p else '0',
            "%s" % (",".join(sorted(three_atype))) if (junction[0], junction[2]) in three_p else '0',
            '\t'.join(additional_fields)])
        snaptron_id+=1
