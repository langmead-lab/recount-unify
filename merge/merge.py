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
MOTIF_COL=6
SAMPLE_COL=7
SAMPLES_COL=8
#separate file columns
FILE_SAMPLE_ID_COL=1

#get the next line, either from the file buffer or a read a new line from the next filehandle
#if there is no more, than that file is done and we'll read another from the heap higher up the stack
def get_record(args, fhs, filebuf, heap, current_chrm, files, last_col, next_fhs):
    i = next_fhs
    fh = fhs[i]
    fields = []
    filedone = False
    if filebuf[i][0] == '' and fh is not None:
        fields = filebuf[i] = fh.readline().rstrip().split('\t')[:last_col]
        #skip header
        if fields[0] == 'chrom':
            fields = filebuf[i] = fh.readline().rstrip().split('\t')[:last_col]
        if filebuf[i][0] == '': fhs[i] = None
    elif filebuf[i][0] != '':
        fields = filebuf[i]
    if len(fields) > 0 and fields[0] != '':
        current_chrm = fields[CHRM_COL]
        if not args.append_samples:
            fields.append(files[i][FILE_SAMPLE_ID_COL])
        heapq.heappush(heap, (fields[CHRM_COL], int(fields[START_COL]), int(fields[END_COL]), (fields, next_fhs)))
        #since we pushed one on the heap we can read another
        if fh is not None:
            filebuf[i] = fh.readline().rstrip().split('\t')[:last_col]
            #check to see if we've exhausted the file
            if filebuf[i][0] == '': fhs[i] = None
    return current_chrm

#either read all of the filehandles (first time) or just the next_fhs
def read_from_sources(args, fhs, filebuf, heap, current_chrm, files, last_col, next_fhs=-1):
    #read 1 line from each source file looking for 1) same chromosome and 2) file EOF
    current_chrm = None
    if next_fhs >= 0:
        current_chrm = get_record(args, fhs, filebuf, heap, current_chrm, files, last_col, next_fhs)
    else:
        for (i, fh) in enumerate(fhs):
            current_chrm = get_record(args, fhs, filebuf, heap, current_chrm, files, last_col, i)
    return current_chrm


def merge(args):
    offset = 0
    scol = SAMPLES_COL
    if args.append_samples:
        offset = 2
        scol = SAMPLE_COL
    motif_col = MOTIF_COL - offset
    last_col = MOTIF_COL + 1
    samples_col = scol - offset
    strand_col = STRAND_COL - offset
    files = []
    with open(args.list_file, "rb") as fin:
        files = [f.rstrip().split('\t') for f in list(fin)]
    if len(args.existing_sj_db) > 0:
        #this only works in --append-samples mode
        files.append([args.existing_sj_db])
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
    next_fh = -1
    while(True):
        #first time (next_fh == -1) read a whole n worth of lines, otherwise just read one more line from the last source that was popped
        current_chrm = read_from_sources(args, fhs, filebuf, heap, current_chrm, files, last_col, next_fh)
        #no more data, end
        #no more data for this chromosome, possibly no more data period, but do another read to find out
        if len(heap) == 0 or heap[0] == '' or heap[0] is None:
            break
        (current, next_fh) = heapq.heappop(heap)[3]
        #print current
        #format sampleID:coverage for Snaptron
        if not args.append_samples:
            current.append(",%s:%s" % (current[SAMPLE_COL],current[COVERAGE_COL]))
        #check for previous junction in case we have to actually merge a junction from multiple sources
        if len(previous) > 0:
            #different junction, print the previous one
            if current[CHRM_COL] != previous[CHRM_COL] or current[START_COL] != previous[START_COL] or current[END_COL] != previous[END_COL]:
                p = previous[:END_COL+1]
                p.extend([previous[strand_col],previous[motif_col],previous[samples_col]])
                sys.stdout.write("%s\n" % ("\t".join(p)))
                previous = current
            #same junction add to the sample IDs/coverages
            else:
                previous[samples_col] += current[samples_col]
        else:
            previous = current
        #sys.stderr.write("FILBUF_SIZE\t"+str(sys.getsizeof(filebuf))+'\t')
        #sys.stderr.write("HEAP_SIZE\t"+str(sys.getsizeof(heap))+'\n')
    #last one
    #sys.stderr.write("FILBUF_SIZE\t"+str(sys.getsizeof(filebuf))+'\t')
    #sys.stderr.write("HEAP_SIZE\t"+str(sys.getsizeof(heap))+'\n')
    if len(previous) > 0:
        p = previous[:END_COL+1]
        p.extend([previous[strand_col],previous[motif_col],previous[samples_col]])
        sys.stdout.write("%s\n" % ("\t".join(p)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='junction merge')
    parser.add_argument('--list-file', metavar='/path/to/manifest_file_of_per_sample_sj_files', type=str, default=None, help='path to a file with list of junction files to aggregate and optionally their sample IDs (rail_ids) and study names)')
    parser.add_argument('--existing-sj-db', metavar='/path/to/existing_db_of_junctions', type=str, default="", help='allows to merge with an existing database of junctionss (assumes same format as output from --append-samples mode), default is empty string')
    parser.add_argument('--gzip', action='store_const', const=True, default=False, help='input files are gzipped')
    parser.add_argument('--append-samples', action='store_const', const=True, default=False, help='set this if aggregating beyond sample-level results')
    args = parser.parse_args()
    
    if len(args.existing_sj_db) > 0 and not args.append_samples:
        sys.stderr.write("ERROR: --existing-sj-db option only works when --append-samples is also specified, exiting!\n")
        sys.exit(-1)

    merge(args)
    
