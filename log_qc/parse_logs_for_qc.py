#!/usr/bin/env python3.6
import sys
import os
import glob
import re

STAR_LOG_SUFFIX = 'align'
KALLISTO_LOG_SUFFIX = 'kallisto'
FC_LOG_SUFFIXES = ['gene_fc_count_all','gene_fc_count_unique','exon_fc_count_all','exon_fc_count_unique']
BAMCOUNT_AUC_SUFFIXES = ['bamcount_auc','bamcount_frag']
SEQTK_SUFFIX='fastq_check'

STAR_PATTERN = re.compile(r'^\s*([\w\s]+)\s\|\t(\d)(\.\d)?')
KALLISTO_PATTERNS = [re.compile(r'[quant]\s+processed\s+([\d,]+)\s+reads,\s+([\d,]+)\s+reads\s+pseudoaligned'), 
FC_PATTERNS = [re.compile(r'(Total\s+[^\s]+)\s*:\s*(\d+)\s'), re.compile(r'(Successfully\s+assigned\s+[^\s]+)\s*:\s*([^\s]+)')]
BAMCOUNT_PATTERN = re.compile(r'^([\t]+)\t(\d+)$')
SEQTK_PATTERN = re.compile(r'^ALL\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)\s+(\d+\.\d+)')

SUFFIXES = {
    'star':[[STAR_LOG_SUFFIX], [STAR_PATTERN]], 
    'kallisto':[[KALLISTO_LOG_SUFFIX], KALLISTO_PATTERNS], 
    'bamcount':[BAMCOUNT_AUC_SUFFIXES, [BAMCOUNT_PATTERN]], 
    'featurecounts':[FC_LOG_SUFFIXES, FC_PATTERNS], 
    'seqtk':[[SEQTK_SUFFIX],[SEQTK_PATTERN]] }

LOGS=['kallisto','star','featurecounts']
TSVS=['bamcount','seqtk']

log_files = glob.glob("%s/**/*.log" % (top_dir), recursive=True)
tsv_files = glob.glob("%s/**/*.tsv*" % (top_dir), recursive=True)




#cat sample_name.auc.tsv for:
#ALL_READS_ANNOTATED_BASES       8497916829
#UNIQUE_READS_ANNOTATED_BASES    7597186502
#ALL_READS_ALL_BASES     9160708402
#UNIQUE_READS_ALL_BASES  8118038679

#parse FC's output for total number, number assigned
#parse seqTK's output


