#!/usr/bin/env python3.6
import sys
import os
import glob
import re
from operator import itemgetter
import subprocess

FILE_FIELD_SEP = '.'
FILE_PREFIX_SEP = '_'
FILE_PREFIX_FIELD_IDX = 0

STAR_SUFFIX = 'align.log'
KALLISTO_SUFFIX = 'kallisto.log'
FC_SUFFIXES = ['gene_fc_count_all.log','gene_fc_count_unique.log','exon_fc_count_all.log','exon_fc_count_unique.log']
SPLIT_LINE_SUFFIXES_MAP = set(FC_SUFFIXES)
BAMCOUNT_SUFFIXES = ['bamcount_auc.tsv','bamcount_frag.tsv']
SEQTK_SUFFIX='fastq_check.tsv.zst'
SPLIT_LINE_SUFFIXES_MAP.add(SEQTK_SUFFIX)
SPLIT_LINE_SUFFIXES_MAP.add(KALLISTO_SUFFIX)

STAR_PATTERN = re.compile(r'^\s*([^\|]+)\s+\|\t(\d+(\.\d+)?)')
KALLISTO_PATTERN = re.compile(r'\[quant\]\s+(processed)\s+([\d,]+)\s+reads,\s+([\d,]+)\s+reads\s+(pseudoaligned)\n\[quant\]\s+(estimated\s+average\s+fragment\s+length):\s+(\d+(\.\d+)?)')
#featurecounts has only 4 fields we want, but spread across 2 lines (consecutively)
FC_PATTERN = re.compile(r'(Total)\s+[^\s]+\s*:\s*(\d+).+\n.+Successfully\s+(assigned)\s+[^\s]+\s*:\s*([^\s]+)')
BAMCOUNT_AUC_PATTERN = re.compile(r'^([^\t]+)\t(\d+)$')
BAMCOUNT_FRAG_PATTERN = re.compile(r'^STAT\t([^\t]+)\t(\d+(\.\d+)?)$')
SEQTK_PATTERN = re.compile(r'ALL\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+).+\nALL\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)')

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
two_group = set([STAR_SUFFIX, BAMCOUNT_SUFFIXES[1], BAMCOUNT_SUFFIXES[0]])
#program patterns that have a 2 label/value groupings (4 group)
four_group = set([KALLISTO_SUFFIX, FC_SUFFIXES[0], FC_SUFFIXES[2], FC_SUFFIXES[2], FC_SUFFIXES[3], SEQTK_SUFFIX])

single_match = four_group

def run_command(cmd_args, cmd_name):
    #dont use shlex.quote, screws up the commandline
    cmd_args = ' '.join(cmd_args)
    sys.stderr.write(cmd_args+'\n')
    try:
        cp = subprocess.run(args=cmd_args, shell=True, check=True, universal_newlines=True) 
    except subprocess.CalledProcessError as cpe:
        sys.stderr.write("error in run_command for command: %s\n" % cmd_args)
        raise cpe

def process_line(line, pattern, suffix, qc):
    match = pattern.search(line)
    matched = False
    if match:
        matched = True
        if suffix in two_group:
            label = match.group(1)
            value = match.group(2)
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
            #we also need a 3rd value (est. avg. frag. len.)
            if suffix == KALLISTO_SUFFIX:
                temp_ = label
                label = value
                value = temp_
                value = value.replace(',','')
                qc[sample]["%s.%s" % (names[suffix],label)] = value
                label = match.group(5)
                value = match.group(6)
            qc[sample]["%s.%s" % (names[suffix],label)] = value
    return matched


top_dir = sys.argv[1]

log_files = glob.glob("%s/**/*.log" % (top_dir), recursive=True)
log_files.extend(glob.glob("%s/**/*.tsv*" % (top_dir), recursive=True))

qc = {}
#separate map for seqTK since it can have a variable # of fields
stk = {}

sample2study = {}
num_mates = 0

for f in log_files:
    (path, file_) = os.path.split(f)
    fields = file_.split(FILE_FIELD_SEP)
    suffix = FILE_FIELD_SEP.join(fields[1:])
    (sample, study, ref) = fields[FILE_PREFIX_FIELD_IDX].split(FILE_PREFIX_SEP)
    if sample not in qc:
        qc[sample] = {}
        stk[sample] = {}
        sample2study[sample] = study
    if suffix not in patterns:
        continue
    pattern = patterns[suffix]
    #need to decode zstd for seqtk log
    if suffix == SEQTK_SUFFIX:
        f1 = "%s.unc" % file_
        run_command(['zstd','-dc',f,' | egrep -e "^ALL" | cut -f 8,9 > %s' % f1], "zstd_seqtk")
        f = f1
    line = None
    with open(f, "r") as fin:
        line = fin.read()
    #handle seqTK differently given it's special format
    if suffix == SEQTK_SUFFIX:
        #os.unlink(f)
        mate_idx = 0
        lines = line.split('\n')
        for line in lines:
            if len(line) == 0:
                continue
            fields_ = line.rstrip().split('\t')
            stk[sample]['%s.P%d.avgQ' % (names[suffix], mate_idx)] = fields_[0] 
            stk[sample]['%s.P%d.errQ' % (names[suffix], mate_idx)] = fields_[1]
            mate_idx += 1
        num_mates = mate_idx
        continue
    if suffix not in SPLIT_LINE_SUFFIXES_MAP:
        lines = line.split('\n')
    else:
        lines = [line]
    for line in lines:
        matched = process_line(line, pattern, suffix, qc)
        #performance, break if we matched and we're only looking for one
        if matched and suffix in single_match:
            break

header = None
header_keys = []
for sample in qc:
    study = sample2study[sample]
    values = qc[sample]
    stks = stk[sample]
    [sys.stderr.write(x+'\t'+values[x]+'\n') for x in values.keys()]
    if header is None:
        header_keys = sorted(values.keys())
        header = '\t'.join([x.lower() for x in header_keys])
        header += '\t'+'\t'.join(['seqtk.P%d.avgQ\tseqtk.P%d.errQ' % (i,i) for i in range(0,num_mates)])
        sys.stdout.write('study\tsample\t'+header+'\n')
    output = '\t'.join([values[x] for x in header_keys])
    for i in range(0,num_mates):
        v = "NA"
        v2 = "NA"
        k = 'seqtk.P%d.avgQ' % i
        k2 = 'seqtk.P%d.errQ' % i
        #if one key is in, then the other will be as well
        if k in stks:
            v = stks[k]
            v2 = stks[k2]
        output += '\t'+v+'\t'+v2
    sys.stdout.write(study+'\t'+sample+'\t'+output+'\n')

