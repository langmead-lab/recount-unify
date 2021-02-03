#!/usr/bin/env python3.6
#join between 2 files using a single field to join on
#default behavior is to do a "outer join" on both sides (left & right)
#i.e. if a key doesn't exist in the other file, print the existing line with blanks for the other file's fields
#this will also do the cross-product join by default
#(every row in each file which matches every row in the other file will get it's own line in the output)
import sys
import os
import glob
import gzip
import argparse
from collections import Counter
from operator import itemgetter
import re
import csv

parser = argparse.ArgumentParser(description='Parse & join non-SRA source metadata for recount3')
parser.add_argument('--smaller', metavar='/path/to/smaller_of_the_two_files.csv', type=str, help='path to file with lesser data (will be hashed on --join-fieldname field), required.  This will also be the file whose fields appear first in a concatenation')
parser.add_argument('--larger', metavar='/path/to/smaller_of_the_two_files.csv', type=str, help='path to file with lesser data (will be streamed and joined to smaller file on --join-fieldname field), required')
parser.add_argument('--smaller-sep', metavar='separator', type=str, default='\t', help='primary field separator for smaller file [default: "\t"]')
parser.add_argument('--larger-sep', metavar='separator', type=str, default='\t', help='primary field separator for larger file [default: "\t"]')
parser.add_argument('--output-sep', metavar='separator', type=str, default='\t', help='primary field sepaarator for output [default: "\t"]')
parser.add_argument('--join-field', metavar='fieldname_to_join_on', type=str, default=None, help='field needed to join across input files (e.g. the key)')
parser.add_argument('--smaller-matching', action='store_const', const=True, default=False, help='don\'t print out unmated rows from the smaller file')
parser.add_argument('--larger-matching', action='store_const', const=True, default=False, help='don\'t print out unmated rows from the larger file')
args = parser.parse_args()

if args.smaller is None:
    sys.stderr.write("ERROR: no file passed in with --smaller, terminating\n")
    exit(-1)

if args.larger is None:
    sys.stderr.write("ERROR: no file passed in with --larger, terminating\n")
    exit(-1)
    
if args.join_field is None:
    sys.stderr.write("ERROR: argument to --join-field , terminating\n")
    exit(-1)

key = args.join_field
smap = {}
smaller_fields = None
with open(args.smaller, newline='') as delimited_file:
    reader = csv.DictReader(delimited_file, restkey='EXTRAFIELDS', delimiter=args.smaller_sep)
    smaller_fields = reader.fieldnames
    for row in reader:
        if row[key] not in smap:
            smap[row[key]] = []
        row_ = args.output_sep.join(row.values()) 
        smap[row[key]].append(row_)

num_smaller_fields = len(smaller_fields)
smaller_blanks_row = args.output_sep.join([""]*num_smaller_fields)

seen = set()
larger_fields = None
with open(args.larger, newline='') as delimited_file:
    reader = csv.DictReader(delimited_file, restkey='EXTRAFIELDS', delimiter=args.larger_sep)
    larger_fields = reader.fieldnames
    #write out joined header
    sys.stdout.write(args.output_sep.join(smaller_fields) + args.output_sep + args.output_sep.join(larger_fields) + '\n')
    for row in reader:
        row_ = args.output_sep.join(row.values()) 
        if row[key] not in smap:
            if args.larger_matching:
                sys.stderr.write("WARNING: row in larger %s not in smaller %s: %s; skipping it\n" % (args.larger, args.smaller, row_))
            else:
                sys.stderr.write("WARNING: row in larger %s not in smaller %s: %s; padding it\n" % (args.larger, args.smaller, row_))
                sys.stdout.write(smaller_blanks_row + args.output_sep + row_ + '\n')
            continue
        smap_rows = smap[row[key]]
        for smap_row in smap_rows:
            sys.stdout.write(smap_row + args.output_sep + row_ + '\n')
        seen.add(row[key])


#now output all rows in smaller file which didn't match anything in larger
num_larger_fields = len(larger_fields)
larger_blanks_row = args.output_sep.join([""]*num_larger_fields)
for k in smap.keys():
    if k in seen:
        continue
    smap_rows = smap[k]
    for smap_row in smap_rows:
        if args.smaller_matching:
            sys.stderr.write("WARNING: row in smaller %s not in larger %s: %s; skipping it\n" % (args.smaller, args.larger, smap_row))
        else:
            sys.stderr.write("WARNING: row in smaller %s not in larger %s: %s; padding it\n" % (args.smaller, args.larger, smap_row))
            sys.stdout.write(smap_row + args.output_sep + larger_blanks_row + '\n')
