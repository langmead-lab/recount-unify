#!/usr/bin/env python3.6
import sys
import os
import glob
import re
from operator import itemgetter

FILE_FIELD_SEP = '.'
FILE_PREFIX_SEP = '_'
FILE_PREFIX_FIELD_IDX = 0

STAR_SUFFIX = 'align.log'
KALLISTO_SUFFIX = 'kallisto.log'
FC_SUFFIXES = ['gene_fc_count_all.log','gene_fc_count_unique.log','exon_fc_count_all.log','exon_fc_count_unique.log']
FC_SUFFIXES_MAP = set(FC_SUFFIXES)
BAMCOUNT_SUFFIXES = ['bamcount_auc.tsv','bamcount_frag.tsv']
SEQTK_SUFFIX='fastq_check.tsv.zst'

STAR_PATTERN = re.compile(r'^\s*([^\|]+)\s+\|\t(\d+(\.\d+)?)')
KALLISTO_PATTERN = re.compile(r'\[quant\]\s+(processed)\s+([\d,]+)\s+reads,\s+([\d,]+)\s+reads\s+(pseudoaligned)')
#featurecounts has only 4 fields we want, but spread across 2 lines (consecutively)
#FC_PATTERN = re.compile(r'(Total\s+[^\s]+)\s*:\s*(\d+).+\n.+(Successfully\s+assigned)\s+[^\s]+\s*:\s*([^\s]+)')
FC_PATTERN = re.compile(r'(Total)\s+[^\s]+\s*:\s*(\d+).+\n.+Successfully\s+(assigned)\s+[^\s]+\s*:\s*([^\s]+)')
BAMCOUNT_AUC_PATTERN = re.compile(r'^([^\t]+)\t(\d+)$')
BAMCOUNT_FRAG_PATTERN = re.compile(r'^STAT\t([^\t]+)\t(\d+(\.\d+)?)$')
SEQTK_PATTERN = re.compile(r'^ALL\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)')

patterns = {
    STAR_SUFFIX:STAR_PATTERN, 
    KALLISTO_SUFFIX:KALLISTO_PATTERN, 
    BAMCOUNT_SUFFIXES[0]:BAMCOUNT_AUC_PATTERN, 
    BAMCOUNT_SUFFIXES[1]:BAMCOUNT_FRAG_PATTERN, 
    FC_SUFFIXES[0]:FC_PATTERN,
    FC_SUFFIXES[1]:FC_PATTERN,
    FC_SUFFIXES[2]:FC_PATTERN,
    FC_SUFFIXES[3]:FC_PATTERN,
    SEQTK_SUFFIX:SEQTK_PATTERN }

#get nice names for each program, keyed by their suffix
names = {s:s.split(FILE_FIELD_SEP)[0] for s in patterns.keys()}
#a little fixup for bamcount and seqtk naming
names[BAMCOUNT_SUFFIXES[0]]='bc_auc'
names[BAMCOUNT_SUFFIXES[1]]='bc_frag'
names[SEQTK_SUFFIX]='seqtk'
names[STAR_SUFFIX]='star'

#program patterns that have a label/value (2 group)
two_group = set([STAR_SUFFIX, BAMCOUNT_SUFFIXES[1], BAMCOUNT_SUFFIXES[0], SEQTK_SUFFIX])
#program patterns that have a 2 label/value groupings (4 group)
four_group = set([KALLISTO_SUFFIX, FC_SUFFIXES[0], FC_SUFFIXES[2], FC_SUFFIXES[2], FC_SUFFIXES[3]])

top_dir = sys.argv[1]

log_files = glob.glob("%s/**/*.log" % (top_dir), recursive=True)
log_files.extend(glob.glob("%s/**/*.tsv*" % (top_dir), recursive=True))

qc = {}

def process_line(line, pattern, suffix, qc):
    match = pattern.search(line)
    if match:
        if suffix in two_group:
            label = match.group(1)
            value = match.group(2)
            if suffix == SEQTK_SUFFIX:
                value_ = label
                label = 'avgQ'
                #use the nice program name.field_label for key
                qc[sample]["%s.%s" % (names[suffix],label)] = value_
                label = 'errQ'
                qc[sample]["%s.%s" % (names[suffix],label)] = value
            else:
                qc[sample]["%s.%s" % (names[suffix],label)] = value
        else:
            label = match.group(1)
            value = match.group(2)
            if suffix == KALLISTO_SUFFIX:
                value = value.replace(',','')
            qc[sample]["%s.%s" % (names[suffix],label)] = value
            label = match.group(3)
            value = match.group(4)
            #Kallisto's pseudoaligned label vs. value is swapped
            if suffix == KALLISTO_SUFFIX:
                temp_ = label
                label = value
                value = temp_
                value = value.replace(',','')
            qc[sample]["%s.%s" % (names[suffix],label)] = value

for f in log_files:
    (path, file_) = os.path.split(f)
    fields = file_.split(FILE_FIELD_SEP)
    suffix = FILE_FIELD_SEP.join(fields[1:])
    (sample, study, ref) = fields[FILE_PREFIX_FIELD_IDX].split(FILE_PREFIX_SEP)
    if sample not in qc:
        qc[sample] = {}
    if suffix not in patterns:
        continue
    pattern = patterns[suffix]
    #TODO need to decode zstd for seqtk log
    if suffix == SEQTK_SUFFIX:
        continue
    line = None
    with open(f, "r") as fin:
        line = fin.read()
    if suffix not in FC_SUFFIXES_MAP:
        lines = line.split('\n')
    else:
        lines = [line]
    for line in lines:
        process_line(line, pattern, suffix, qc)

header = None
for sample in qc:
    values = qc[sample]
    [sys.stderr.write(x+'\t'+values[x]+'\n') for x in values.keys()]
    if header is None:
        header = '\t'.join(sorted([x.lower() for x in values.keys()]))
        sys.stdout.write(header+'\n')
    output = '\t'.join([values[x] for x in sorted(values.keys())])
    sys.stdout.write(output+'\n')

