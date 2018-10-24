#!/bin/env python
import sys
import heapq

CHRM_COL=0
START_COL=1
END_COL=2
COVERAGE_COL=4
STRAND_COL=5
SAMPLE_COL=6
SAMPLES_COL=7
LAST_COL=STRAND_COL+1
#separate file columns
FILE_SAMPLE_ID_COL=1


def read_from_sources(fhs, filebuf, heap, current_chrm, files):
    #read 1 line from each source file looking for 1) same chromosome and 2) file EOF
    for (i,fin) in enumerate(fhs):
        fields = []
        filedone = False
        if filebuf[i][0] == '' and fin is not None:
            fields = filebuf[i] = fin.readline().rstrip().split('\t')[:LAST_COL]
            #skip header
            if fields[0] == 'chrom':
                fields = filebuf[i] = fin.readline().rstrip().split('\t')[:LAST_COL]
            filedone = filebuf[i][0] == ''
        elif filebuf[i][0] != '':
            fields = filebuf[i]
        if len(fields) > 0 and (current_chrm is None or fields[CHRM_COL] == current_chrm):
            current_chrm = fields[CHRM_COL]
            #make sure we trak the sample ID
            fields.append(files[i][FILE_SAMPLE_ID_COL])
            heapq.heappush(heap, (int(fields[START_COL]), fields))
            if fin is not None:
                filebuf[i] = fin.readline().rstrip().split('\t')[:LAST_COL]
                filedone = filebuf[i][0] == ''
        if filedone: fhs[i] = None
    #if this is None, means we're at the end of all files
    return current_chrm


def merge(sample_file_list):
    files = []
    with open(sample_file_list, "rb") as fin:
        files = [f.rstrip().split('\t') for f in list(fin)]
    fhs = [open(f[0],"rb") for f in files]
    n = len(fhs)
    filebuf = [['']]*n
    heap = []
    current_chrm = None
    previous = []
    while(True):
        current_chrm = read_from_sources(fhs, filebuf, heap, current_chrm, files)
        #no more data, end
        if current_chrm is None: break
        #no more data for this chromosome, possibly no more data period, but do another read to find out
        if len(heap) == 0 or heap[0] is None:
            current_chrm = None
            continue
        current = heapq.heappop(heap)[1]
        #print current
        #format sampleID:coverage for Snaptron
        current.append(",%s:%s" % (current[SAMPLE_COL],current[COVERAGE_COL]))
        #check for previous junction in case we have to actually merge a junction from multiple sources
        if len(previous) > 0:
            #different junction, print the previous one
            if current[CHRM_COL] != previous[CHRM_COL] or current[START_COL] != previous[START_COL] or current[END_COL] != previous[END_COL]:
                p = previous[:END_COL+1]
                p.extend([previous[STRAND_COL],previous[SAMPLES_COL]])
                sys.stdout.write("%s\n" % ("\t".join(p)))
                previous = current
            #same junction add to the sample IDs/coverages
            else:
                previous[SAMPLES_COL] += current[SAMPLES_COL]
        else:
            previous = current
    #last one
    if len(previous) > 0:
        p = previous[:END_COL+1]
        p.extend([previous[STRAND_COL],previous[SAMPLES_COL]])
        sys.stdout.write("%s\n" % ("\t".join(p)))

if __name__ == '__main__':
    sample_file_list = sys.argv[1]
    merge(sample_file_list)
    
