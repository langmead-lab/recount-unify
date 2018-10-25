#!/bin/env python
import sys
import heapq
import argparse
import gzip

CHRM_COL=0
START_COL=1
END_COL=2
COVERAGE_COL=4
STRAND_COL=5
SAMPLE_COL=6
SAMPLES_COL=7
#separate file columns
FILE_SAMPLE_ID_COL=1


def read_from_sources(fhs, filebuf, heap, current_chrm, files, last_col):
    #read 1 line from each source file looking for 1) same chromosome and 2) file EOF
    for (i,fin) in enumerate(fhs):
        fields = []
        filedone = False
        if filebuf[i][0] == '' and fin is not None:
            fields = filebuf[i] = fin.readline().rstrip().split('\t')[:last_col]
            #skip header
            if fields[0] == 'chrom':
                fields = filebuf[i] = fin.readline().rstrip().split('\t')[:last_col]
            filedone = filebuf[i][0] == ''
        elif filebuf[i][0] != '':
            fields = filebuf[i]
        if len(fields) > 0 and (current_chrm is None or fields[CHRM_COL] == current_chrm):
            current_chrm = fields[CHRM_COL]
            #make sure we trak the sample ID
            fields.append(files[i][FILE_SAMPLE_ID_COL])
            heapq.heappush(heap, (int(fields[START_COL]), fields))
            if fin is not None:
                filebuf[i] = fin.readline().rstrip().split('\t')[:last_col]
                filedone = filebuf[i][0] == ''
        if filedone: fhs[i] = None
    #if this is None, means we're at the end of all files
    return current_chrm


def merge(args):
    last_col = STRAND_COL+1
    samples_col = SAMPLES_COL
    if args.append_samples:
        last_col = STRAND_COL
        samples_col = COVERAGE_COL
    files = []
    with open(args.list_file, "rb") as fin:
        files = [f.rstrip().split('\t') for f in list(fin)]
    fhs = []
    if args.gzip:
        fhs = [gzip.open(f[0],"rb") for f in files]
    else:
        fhs = [open(f[0],"rb") for f in files]
    n = len(fhs)
    filebuf = [['']]*n
    heap = []
    current_chrm = None
    previous = []
    while(True):
        current_chrm = read_from_sources(fhs, filebuf, heap, current_chrm, files, last_col)
        #no more data, end
        if current_chrm is None: break
        #no more data for this chromosome, possibly no more data period, but do another read to find out
        if len(heap) == 0 or heap[0] is None:
            current_chrm = None
            continue
        current = heapq.heappop(heap)[1]
        #print current
        #format sampleID:coverage for Snaptron
        if not args.append_samples:
            current.append(",%s:%s" % (current[SAMPLE_COL],current[COVERAGE_COL]))
        #check for previous junction in case we have to actually merge a junction from multiple sources
        if len(previous) > 0:
            #different junction, print the previous one
            if current[CHRM_COL] != previous[CHRM_COL] or current[START_COL] != previous[START_COL] or current[END_COL] != previous[END_COL]:
                p = previous[:END_COL+1]
                p.extend([previous[STRAND_COL],previous[samples_col]])
                sys.stdout.write("%s\n" % ("\t".join(p)))
                previous = current
            #same junction add to the sample IDs/coverages
            else:
                previous[samples_col] += current[samples_col]
        else:
            previous = current
    #last one
    if len(previous) > 0:
        p = previous[:END_COL+1]
        p.extend([previous[STRAND_COL],previous[samples_col]])
        sys.stdout.write("%s\n" % ("\t".join(p)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='junction merge')
    parser.add_argument('--list-file', metavar='/path/to/file_with_queries', type=str, default=None, help='path to a file with list of junction files to aggregate and optionally their sample IDs (rail_ids) and study names)')
    parser.add_argument('--gzip', action='store_const', const=True, default=False, help='input files are gzipped')
    parser.add_argument('--append-samples', action='store_const', const=True, default=False, help='set this if aggregating beyond sample-level results')
    args = parser.parse_args()
    merge(args)
    
