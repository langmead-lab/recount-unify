#!/usr/bin/env python3.6
import sys
import os
import glob
import re
from operator import itemgetter
import subprocess
import argparse
from collections import Counter
import patroller
import pickle
import gzip

def load_cpickle_file(filepath, compressed=False):
    ds = None
    if os.path.exists(filepath):
        if compressed:
            with gzip.GzipFile(filepath,"rb") as f:
                ds=pickle.load(f)
        else: 
            with open(filepath,"rb") as f:
                ds=pickle.load(f)
    return ds

def store_cpickle_file(filepath, ds, compress=False):
    if not os.path.exists(filepath):
        if compress:
            with gzip.GzipFile(filepath,"wb") as f:
                pickle.dump(ds,f,pickle.HIGHEST_PROTOCOL)
        else:
            with open(filepath,"wb") as f:
                pickle.dump(ds,f,pickle.HIGHEST_PROTOCOL)
        return True
    return False

#always double STAR columns to be consistent for those samples that have the 3rd FASTQ file from SRA
#(which leads to having 2 columns for every STAR aligner stat)
STAR_twice = True

FILE_FIELD_SEP = '.'
#FILE_PREFIX_SEP = '_'
FILE_PREFIX_SEP = '!'
FILE_PREFIX_FIELD_IDX = 0

IDXSTATS_SUFFIX='idxstats'
STAR_SUFFIX = 'align.log'
FC_SUFFIXES = ['gene_fc_count_all.log','gene_fc_count_unique.log','exon_fc_count_all.log','exon_fc_count_unique.log']
SPLIT_LINE_SUFFIXES_MAP = set(FC_SUFFIXES)
BAMCOUNT_SUFFIXES = ['bamcount_auc.tsv','bamcount_frag.tsv']
SEQTK_SUFFIX='fastq_check.tsv.zst'
SPLIT_LINE_SUFFIXES_MAP.add(SEQTK_SUFFIX)

STAR_PATTERN = re.compile(r'^\s*([^\|]+)\s+\|\t(\d+(\.\d+)?)')
#featurecounts has only 4 fields we want, but spread across 2 lines (consecutively)
FC_PATTERN = re.compile(r'(Total)\s+[^\s]+\s*:\s*(\d+).+\n.+Successfully\s+(assigned)\s+[^\s]+\s*:\s*([^\s]+)')
BAMCOUNT_AUC_PATTERN = re.compile(r'^([^\t]+)\t(\d+)$')
BAMCOUNT_FRAG_PATTERN = re.compile(r'^STAT\t([^\t]+)\t(\d+(\.\d+)?)$')
SEQTK_PATTERN = re.compile(r'ALL\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+).+\nALL\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)')
#get mapping rate for chrs: X, Y, and M for QC
IDXSTATS_PATTERN=re.compile(r'^([^\*][^\t]*)\t\d+\t(\d+)\t')
IDXSTATS_CHRS=set(['chrM','chrY','chrX'])

patterns = {
    IDXSTATS_SUFFIX:IDXSTATS_PATTERN,
    STAR_SUFFIX:STAR_PATTERN, 
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
names[IDXSTATS_SUFFIX]='aligned_reads%'

#program patterns that have a label/value (2 group)
two_group = set([STAR_SUFFIX, BAMCOUNT_SUFFIXES[1], BAMCOUNT_SUFFIXES[0], IDXSTATS_SUFFIX])
#program patterns that have a 2 label/value groupings (4 group)
#four_group = set([KALLISTO_SUFFIX, FC_SUFFIXES[0], FC_SUFFIXES[2], FC_SUFFIXES[2], FC_SUFFIXES[3], SEQTK_SUFFIX])
four_group = set([FC_SUFFIXES[0], FC_SUFFIXES[2], FC_SUFFIXES[2], FC_SUFFIXES[3], SEQTK_SUFFIX])

single_match = four_group

#keep track of total set of unique field names as keys
keys = set()

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
            if suffix == IDXSTATS_SUFFIX:
                qc[sample]['idxstats.all_mapped_reads'] += int(value)
                #if not chrs M, X, or Y don't save them specifically
                if label not in IDXSTATS_CHRS:
                    return matched
            if suffix == STAR_SUFFIX:
                label = label.replace(' ','_')
                key = "%s.%s" % (names[suffix],label)
                #if STAR stat name is already present in the hash
                #add it as a 2nd one
                if key in qc[sample]:
                    key = key + "2"
                #go ahead and add a default value
                else:
                    key2 = key + "2"
                    qc[sample][key2]=0; 
                    keys.add(key2)
            else:
                key = "%s.%s" % (names[suffix],label)
            qc[sample][key] = value
            keys.add(key)
        else:
            label = match.group(1)
            value = match.group(2)
            key = "%s.%s" % (names[suffix],label)
            qc[sample][key] = value
            keys.add(key)
            label = match.group(3)
            value = match.group(4)
            key = "%s.%s" % (names[suffix],label)
            qc[sample][key] = value
            keys.add(key)
    return matched

def load_intron_sums(args):
    sample_map = {}
    with open(args.sample_mapping,'r') as fin:
        sample_map = {line.split('\t')[2]:line.split('\t')[1] for line in fin.read().split('\n')[:-1]}
    with open(args.intron_sums,'r') as fin:
        header = fin.readline().rstrip().split('\t')
        return {sample_map[header[i+1]]:intron_sum for (i,intron_sum) in enumerate(fin.readline().rstrip().split('\t')[1:])}

parser = argparse.ArgumentParser(description='Parse log/summary stat files from monorail run')
#e.g. --incoming-dir /home-1/cwilks3@jhu.edu/storage/recount-pump/destination/geuv_sc.20190206
parser.add_argument('--incoming-dir', metavar='/path/to/dir_containing_pump_processed_files', type=str, default=None, help='the path where recount-pump dumps it\'s procssed files.')
parser.add_argument('--sample-mapping', metavar='/path/to/sample_id_mapping_file.tsv', type=str, default=None, help='path to sample ID mapping to accession [optional, use only if joining intron sums')
parser.add_argument('--intron-sums', metavar='/path/to/intron_sums.tsv', type=str, default=None, help='path to intron sums, one value per sample [optional, but requires the setting of --sample-mapping]')
parser.add_argument('--dont-use-patroller', action='store_const', const=True, default=False, help='if user has already filtered the attempts to be the correct set then they can skip using the patroller to find finished monorail runs')
parser.add_argument('--use-pickled', action='store_const', const=True, default=False, help='if re-running on a large set and don\'t want to re-run the find/glob again')
args = parser.parse_args()

intron_sums = {}
if args.intron_sums is not None and args.sample_mapping is None:
    sys.stderr.write("intron sums passed in but no sample mapping, skipping inron sums\n")
elif args.intron_sums is not None:
    intron_sums = load_intron_sums(args) 

intron_sums_len = len(intron_sums)

top_dir = args.incoming_dir

seen = {}
loworders = {}
done = {}
attempts_tracker = {}
runs_by_study_done = Counter()
#only use "seen" here, the rest are just for compatibility
patroller.find_done_runs(args, {}, runs_by_study_done, seen, {}, {})

fp="patroller_find_done.tsv"
#if not os.path.exists(fp):
with open(fp,"w") as fout:
    [fout.write("%s\t%s\t%s\n" % (study, str(fkey), "\t".join([str(x) for x in seen[study][fkey]]))) for study in seen.keys() for fkey in seen[study].keys()]

fp="all_logs_and_tsvs.pkl"
log_files = None
if args.use_pickled:
    log_files = load_cpickle_file(fp)
if log_files is None:
    log_files = glob.glob("%s/**/*.idxstats" % (top_dir), recursive=True)
    log_files.extend(glob.glob("%s/**/*.log" % (top_dir), recursive=True))
    log_files.extend(glob.glob("%s/**/*.tsv*" % (top_dir), recursive=True))
    store_cpickle_file(fp, log_files) 

qc = {}
#separate map for seqTK since it can have a variable # of fields
stk = {}

sample2study = {}
num_mates = 0

dups = set()

fpath_patt = re.compile(r'^((.+)_attempt(\d+))')
fpath_patt2 = re.compile(r'^((.+)_att(\d+))')
for f in log_files:
    (path, file_) = os.path.split(f)
    if file_ in dups:
        continue
    fields = file_.split(FILE_FIELD_SEP)
    if len(fields) < 2:
        sys.stderr.write("%s format not expected, skipping\n" % (f))
        continue
    suffix = FILE_FIELD_SEP.join(fields[1:])
    fields2 = fields[FILE_PREFIX_FIELD_IDX].split(FILE_PREFIX_SEP)
    if len(fields2) < 3:
        sys.stderr.write("%s format not expected, skipping\n" % (f))
        continue
    (sample, study, ref) = fields2[:3]
    m = fpath_patt.search(path)
    if m is None:
        m = fpath_patt2.search(f)
    fkey = m.group(2)
    attempt_num = int(m.group(3))
    #keep track of what study, proj_path, and attempt
    if not args.dont_use_patroller and (study not in seen or fkey not in seen[study] or seen[study][fkey][0] != attempt_num):
        continue
    if sample not in qc:
        qc[sample] = {}
        stk[sample] = {}
        sample2study[sample] = study
        qc[sample]['idxstats.all_mapped_reads'] = 0
        #keys.add('idxstats.all_mapped_reads')
    if suffix not in patterns:
        continue
    dups.add(file_)
    pattern = patterns[suffix]
    #need to decode zstd for seqtk log
    #TODO re-enable SEQTK
    if suffix == SEQTK_SUFFIX:
        continue
        #f1 = "%s.unc" % file_
        #run_command(['zstd','-dc',f,' | egrep -e "^ALL" | cut -f 8,9 > %s' % f1], "zstd_seqtk")
        #f = f1
    line = None
    with open(f, "r") as fin:
        line = fin.read()
    #handle seqTK differently given it's special format
    if suffix == SEQTK_SUFFIX:
        os.unlink(f)
        mate_idx = 0
        lines = line.split('\n')
        for line in lines:
            if len(line) == 0:
                continue
            fields_ = line.rstrip().split('\t')
            key = '%s.P%d.avgQ' % (names[suffix], mate_idx)
            stk[sample][key] = fields_[0] 
            key = '%s.P%d.errQ' % (names[suffix], mate_idx)
            stk[sample][key] = fields_[1]
            #keys.add(key)
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

ratio_cols = [
            ['bc_auc.all_%','bc_auc.ALL_READS_ANNOTATED_BASES','bc_auc.ALL_READS_ALL_BASES'],
            ['bc_auc.unique_%','bc_auc.UNIQUE_READS_ANNOTATED_BASES','bc_auc.UNIQUE_READS_ALL_BASES'],
            ['gene_fc.all_%','gene_fc_count_all.assigned','star.All_mapped_reads'],
            ['gene_fc.unique_%','gene_fc_count_unique.assigned','star.Uniquely_mapped_reads_number'],
            ['exon_fc.all_%','exon_fc_count_all.assigned','star.All_mapped_reads'],
            ['exon_fc.unique_%','exon_fc_count_unique.assigned','star.Uniquely_mapped_reads_number']
]

keys.update([x[0] for x in ratio_cols])

#these might not be present
optional_ratio_cols = ['gene_fc_count_all.assigned','gene_fc_count_unique.assigned','exon_fc_count_all.assigned','exon_fc_count_unique.assigned']

keys.add('star.All_mapped_reads')
if STAR_twice:
    keys.add('star.All_mapped_reads2')
if intron_sums_len > 0:
    keys.add('intron_sum_%')
    keys.add('intron_sum')
header = None
header_keys = sorted(keys)
#header_keys = [col.sub(' ','_') for col in sorted(keys)]
header = '\t'.join([x.lower() for x in header_keys])
#header += '\t'+'\t'.join(['seqtk.P%d.avgQ\tseqtk.P%d.errQ' % (i,i) for i in range(0,num_mates)])
sys.stdout.write('study\tsample\t'+header+'\n')
for sample in qc:
    study = sample2study[sample]
    values = qc[sample]
    stks = stk[sample]
    output = []
    #make sure all the keys are set
    for x in header_keys:
        if x not in values:
            values[x] = '0'
    if 'gene_fc_count_all.assigned' not in values or values['gene_fc_count_all.assigned'] == '0':
        for k in optional_ratio_cols:
            values[k] = '0'
    values['star.All_mapped_reads'] = str(int(values['star.Uniquely_mapped_reads_number']) + \
                                        int(values['star.Number_of_reads_mapped_to_multiple_loci']))
    if STAR_twice:
        values['star.All_mapped_reads2'] = str(int(values['star.Uniquely_mapped_reads_number2']) + \
                                            int(values['star.Number_of_reads_mapped_to_multiple_loci2']))

    #do %'s of idxstats
    denom_ = values['idxstats.all_mapped_reads']
    for chrm in IDXSTATS_CHRS:
        key_ = '%s.%s' % (names[IDXSTATS_SUFFIX], chrm)
        if denom_ > 0 and key_ in values:
            values[key_] = round(100*(float(values[key_])/denom_),2)
        else:
            values[key_] = -1.0

    bc_auc_all_all = 0
    if 'bc_auc.ALL_READS_ALL_BASES' in values:
        bc_auc_all_all = int(values['bc_auc.ALL_READS_ALL_BASES'])
    if intron_sums_len > 0:
        values['intron_sum_%'] = 0
        values['intron_sum'] = 0
        if sample in intron_sums and bc_auc_all_all > 0:
            values['intron_sum_%'] = str(round(100*(float(intron_sums[sample])/bc_auc_all_all),2))
            values['intron_sum'] = intron_sums[sample]

    #dont need this as an output value will just be confusing
    #as it measures the total number of read mappings which can include duplicates of reads and is 2x the # of fragments 
    del(values['idxstats.all_mapped_reads'])
    #count of fragments vs. total input from STAR
    num_frags = values['bc_frag.COUNT']
    if int(num_frags) != int(values['star.All_mapped_reads']):
        sys.stderr.write("WARNING\tNUM_FRAGS != STAR MAPPED READS for %s: %s vs. %s\n" % (sample, str(num_frags),str(values['star.All_mapped_reads']))) 
    total_input_frags = values['star.Number_of_input_reads']
    #values['mapped fragments/total input %'] = str(int(100*round(int(num_frags) / int(total_input_frags),2)))
    values.update({n[0]:'0' for n in ratio_cols})
    values.update({n[0]:str(round(100*int(values[n[1]])/int(values[n[2]]),2)) for n in ratio_cols if int(values[n[2]]) != 0.0})
    #[sys.stderr.write(x+'\t'+str(values[x])+'\n') for x in header_keys]
    output = '\t'.join([str(values[x]) for x in header_keys])
    sys.stdout.write(study+'\t'+sample+'\t'+output+'\n')
    continue
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
    #sys.stdout.write(study+'\t'+sample+'\t'+output+'\n')
