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
parser.add_argument('--smaller-matching', action='store_const', const=True, default=False, help='don\'t print out unmatched rows from the smaller file')
parser.add_argument('--larger-matching', action='store_const', const=True, default=False, help='don\'t print out unmatched rows from the larger file')

#parser.add_argument('--keep-case', action='store_const', const=True, default=False, help='don\'t make all field names lowercase (for R) default: lower all fieldnames')

#parser.add_argument('--print-duplicate-join-fields', action='store_const', const=True, default=False, help='don\'t suppress the printing of additional (duplicate) join fields from all the files, default is to only print one of them (from --smaller file)')

parser.add_argument('--join-field', metavar='fieldname_to_join_on', type=str, help='field needed to join across input files (if the same) or for the smaller file (if different), required')
parser.add_argument('--larger-field', metavar='fieldname_to_join_on', type=str, default=None, help='field needed to join the larger input file on (e.g. the key), defaults to the value for --smaller-key')

args = parser.parse_args()

if args.smaller is None:
    sys.stderr.write("ERROR: no file passed in with --smaller, terminating\n")
    exit(-1)

if args.larger is None:
    sys.stderr.write("ERROR: no file passed in with --larger, terminating\n")
    exit(-1)
    
if args.join_field is None:
    sys.stderr.write("ERROR: argument to --join-field, terminating\n")
    exit(-1)

smaller_key = args.join_field
larger_key = args.larger_field
if larger_key is None:
    larger_key = smaller_key

smap = {}
smaller_fields_deduped = []
smaller_fields_set = set()
smaller_idxs = []
with open(args.smaller, newline='') as delimited_file:
    reader = csv.DictReader(delimited_file, restkey='EXTRAFIELDS', delimiter=args.smaller_sep)
    for (i,field) in enumerate(reader.fieldnames):
        field_purged = field.replace('\n', ' ').replace(' ', '_')
        #TODO make lower an option
        field_purged = field_purged.lower()
        if field_purged not in smaller_fields_set:
            smaller_idxs.append(i)
            smaller_fields_set.add(field_purged)
            smaller_fields_deduped.append(field_purged)
    for row in reader:
        if row[smaller_key] not in smap:
            smap[row[smaller_key]] = []
        #still need original casing for matching reader's row fields
        fields = []
        for i in smaller_idxs:
            field = row[reader.fieldnames[i]]
            if field is not None:
                field = field.replace('\n',' ')
            else:
                field = ""
            fields.append(field)
        row_ = args.output_sep.join(fields)
        #use original key value for matching with larger file
        smap[row[smaller_key]].append(row_)

num_smaller_fields = len(smaller_fields_deduped)
smaller_blanks_row = args.output_sep.join([""]*num_smaller_fields)
    
seen = set()
larger_fields = []
larger_fields_deduped = []
larger_fields_set = set()
larger_idxs = []
with open(args.larger, newline='') as delimited_file:
    reader = csv.DictReader(delimited_file, restkey='EXTRAFIELDS', delimiter=args.larger_sep)
    for (i,field) in enumerate(reader.fieldnames):
        field_purged = field.replace('\n', ' ').replace(' ', '_')
        #TODO make lower an option
        field_purged = field_purged.lower()
        if field_purged not in larger_fields_set and field_purged not in smaller_fields_set:
            larger_idxs.append(i)
            larger_fields_deduped.append(field_purged)
            larger_fields_set.add(field_purged)
    #write out joined header
    sys.stdout.write(args.output_sep.join(smaller_fields_deduped) + args.output_sep + args.output_sep.join(larger_fields_deduped) + '\n')
    for row in reader:
        fields = []
        for i in larger_idxs:
            field = row[reader.fieldnames[i]]
            if field is not None:
                field = field.replace('\n',' ')
            else:
                field = ""
            fields.append(field)
        row_ = args.output_sep.join(fields)
        #use original key value for matching with smaller file
        if row[larger_key] not in smap:
            if args.larger_matching:
                sys.stderr.write("WARNING: row in larger %s not in smaller %s: %s; skipping it\n" % (args.larger, args.smaller, row_))
            else:
                sys.stderr.write("WARNING: row in larger %s not in smaller %s: %s; padding it\n" % (args.larger, args.smaller, row_))
                sys.stdout.write(smaller_blanks_row + args.output_sep + row_ + '\n')
            continue
        smap_rows = smap[row[larger_key]]
        for smap_row in smap_rows:
            sys.stdout.write(smap_row + args.output_sep + row_ + '\n')
        seen.add(row[larger_key])
    
#now output all rows in smaller file which didn't match anything in larger
num_larger_fields = len(larger_fields_deduped)
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
