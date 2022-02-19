#!/usr/bin/env python
#input: the gene level output of the rejoin program
#output: produces the gene level sums for each original annotation in separate files

import sys
import re

#mapping file between gene extend coords and gene name/annotation source
mapf = sys.argv[1]
#gene or exon
annot_type = sys.argv[2]
counts_headerF = sys.argv[3]
annot_source = sys.argv[4]
id_mappingF = sys.argv[5]
#e.g. G026, write out a different header for Snaptron which 
#uses the rail_ids instead of the sample accessions/barcodes
snaptron_annot_src = sys.argv[6]
study = None
if len(sys.argv) >= 8:
    study = sys.argv[7]

annot_sources = set(annot_source.split(','))

counts_header = ""
with open(counts_headerF,"r") as fin:
    counts_header = fin.read().rstrip()

#get accessions (or barcodes) as counts header rather than rail_ids
if study is not None:
    with open(snaptron_annot_src+"."+annot_type+"_sums.snaptron_header.%s.tsv" % study,"w") as fout:
        fout.write("gene_id\tbp_length\tchromosome\tstart\tend\t%s\n" % (counts_header))
else:
    with open(snaptron_annot_src+"."+annot_type+"_sums.snaptron_header.tsv","w") as fout:
        fout.write("gene_id\tbp_length\tchromosome\tstart\tend\t%s\n" % (counts_header))

id_mapping = {}
with open(id_mappingF,"r") as fin:
    for line in fin:
        (study_, run, rail_id) = line.rstrip().split('\t')
        id_mapping[rail_id] = run

counts_header_fields = counts_header.split('\t')
#may be using accessions/original IDs already
try: 
    int(counts_header_fields[0])
    counts_header = '\t'.join([id_mapping[rid] for rid in counts_header_fields])
except ValueError:
    pass

annotation_fhs = {}
annot_map = {}

additional_annots_patt = re.compile(r'^((SIRV)|(ERCC))')
def determine_annotation_source(gname, annot_fhs=None):
    gsplit = gname.split('.')
    annot_src = None
    num_splits = len(gsplit)
    found_idx = -1
    for i,potential_src in enumerate(gsplit):
        if potential_src in annot_sources:
            annot_src = potential_src
            del gsplit[i]
            break
    if annot_src is None:
        m = additional_annots_patt.search(gname)
        if m is not None:
            annot_src = m.group(1)
    gene_id = '.'.join(gsplit)
    if annot_src is not None and annot_fhs is not None and annot_src not in annot_fhs:
        if study is not None:
            annot_fhs[annot_src] = open(annot_src+"."+annot_type+"_sums.%s.tsv" % study,"w")
        else:
            annot_fhs[annot_src] = open(annot_src+"."+annot_type+"_sums.tsv","w")
        annot_fhs[annot_src].write("gene_id\tbp_length\tchromosome\tstart\tend\t%s\n" % (counts_header))
    return (gene_id, annot_src)


#check if we're doing mouse
mouse = snaptron_annot_src[0] == 'M'

with open(mapf,"r") as fin:
    for (idx, line) in enumerate(fin):
        line = line.rstrip()
        (chrm, gstart, gend, gname, typenum, strand)  = line.split('\t')
        gstart = str(int(gstart)+1)

        (gene_id, annot_src) = determine_annotation_source(gname, annotation_fhs)
        if annot_src is None:
            continue

        key = "\t".join([chrm, gstart, gend, strand])
        if key not in annot_map:
            annot_map[key] = {}
        if gene_id not in annot_map[key]:
            annot_map[key][gene_id] = []
        annot_map[key][gene_id].append(annot_src)

#avoid writing duplicates
seen = set()
for line in sys.stdin:
    line = line.rstrip()
    fields = line.split('\t')
    (gnames, chrm, gstart, gend, glength, strand) = fields[:6]
    #account for 0-based start coordinate (depending on organism)
    gstart = str(int(gstart)+1)
    counts = "\t".join(fields[6:])
    #output is: Gene   bp_length       chromosome      start   end count1 count2 count3....
    #outstr = "\t".join([str((int(gend) - int(gstart)) + 1), chrm, gstart, gend])+"\t"+counts+"\n"
    outstr = "\t".join([glength, chrm, gstart, gend])+"\t"+counts+"\n"
    key = "\t".join([chrm, gstart, gend, strand])
    
    #now figure out the mapping to all the of the genes and annotations
    genes_and_annotations = []
    #mouse has only one real annoation (M023 for now) 
    if mouse and gnames[-4:] != 'SIRV' and gnames[-4:] != 'ERCC':
        for gname in gnames.split(';'):
            genes_and_annotations.append([gname, snaptron_annot_src])
    else:
        for gname in gnames.split(';'): 
            (gene_id, annot_src) = determine_annotation_source(gname)
            #print "annot_src determination "+str(gname)+":"+str(annot_src)
            if annot_src is None:
                continue
            genes_and_annotations.append([gene_id, annot_src])
    for (gene_id, annot_src) in genes_and_annotations:
        for annot_src in annot_map[key][gene_id]:
            key2 = gene_id + "_" + annot_src
            if key2 not in seen:
                annotation_fhs[annot_src].write(gene_id + "\t" + outstr)
                seen.add(key2)

for fh in annotation_fhs.values():
    fh.close()
