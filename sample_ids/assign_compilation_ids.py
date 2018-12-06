#!/bin/env python2
import sys
import shutil
import argparse
from collections import Counter
import buildIDs

COMPILATION_CODE_BIT_LENGTH=8

def main(args):
    #currently only converts a SRA run accession (e.g. SRR1234567) or a TCGA GDC File ID (UUID) (e.g. bef8b87c-4473-441a-a840-0bd8f2f1bc9e)
    counter = -1
    compilation_code = ""
    if args.compilation_code != -1:
        compilation_code = str(bin(args.compilation_code)[2:].zfill(COMPILATION_CODE_BIT_LENGTH))
    prefix_counters = Counter()
    if args.group_counters_file is not None:
        with open(args.group_counters_file, "r") as fin:
            for line in fin:
                (group, count) = line.rstrip().split('\t')
                prefix_counters[str(group)] = int(count)
    fin = open(args.accessions_file, "rb")
    for line in fin:
        if counter == -1 and not args.no_header:
            counter+=1
            continue
        fields = line.rstrip().split('\t')
        if len(fields) <= args.acc_col:
            sys.stderr.write("ERROR: --acc-col of %d is larger than number of fields in %s, terminating!\n" % (args.acc_col, args.accessions_file))
            sys.exit(-1)
        acc = fields[args.acc_col]
        prefix_str = buildIDs.encodeString(acc, buildIDs.alphaParseSuffix)
        #sys.stderr.write(prefix_str+'\n')
        #use prefix-relative counters
        counter_str = str(bin(prefix_counters[prefix_str])[2:])
        prefix_counters[prefix_str] += 1
        prefix_ctr_str = compilation_code+prefix_str+counter_str
        #sys.stdout.write(str(int(prefix_ctr_str,2))+"\t"+acc+\t"+"\n")
        sys.stdout.write(fields[0]+"\t"+acc+"\t"+str(int(prefix_ctr_str,2))+"\n")
        counter+=1
    shutil.move(args.group_counters_file, "%s.old" % (args.group_counters_file))
    with open(args.group_counters_file, "w") as fout:
        [fout.write("%s\t%s\n" % (x, str(prefix_counters[x]))) for x in prefix_counters.keys()]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assigns IDs for a specified compilation')

    parser.add_argument('--accessions-file', metavar='/path/to/file_with_accessions', type=str, default=None, help='path to file with list of accessions to assign IDs to')
    parser.add_argument('--group-counters-file', metavar='/path/to/file_with_current_group_id_counters', type=str, default=None, help='path to file with the list of groups (low order digits) and their last used id counts')
    parser.add_argument('--acc-col', metavar='integer column #', type=int, default=0, help='column number in --accessions-file where accession is (defaults to 0)')
    parser.add_argument('--compilation-code', metavar='8-bit integer code per compilation', type=int, default=-1, help='used as a prefix to the sample ID to differentiate samples from different compilations when results are merged across compilations; 15 (1111) is reserved for local (non-public) compilations.  Default is not to include the compilation ID in the ID generation.')
    parser.add_argument('--no-header', action='store_const', const=True, default=False, help='if --accessions-file doesn\'t have a header')
    
    args = parser.parse_args()

    main(args)
