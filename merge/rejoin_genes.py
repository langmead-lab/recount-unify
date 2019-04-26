#!/usr/bin/env python
#input: the gene level output of the rejoin program
#output: produces the gene level sums for each original annotation in separate files

import sys

#mapping file between gene extend coords and gene name/annotation source
mapf = sys.argv[1]
#gene or exon
annot_type = sys.argv[2]
counts_headerF = sys.argv[3]

counts_header = ""
with open(counts_headerF,"r") as fin:
    counts_header = fin.read().rstrip()

annot_fhs = {}
annot_map = {}

with open(mapf,"r") as fin:
    for line in fin:
        line = line.rstrip()
        (chrm, gstart, gend, gname, typenum, strand)  = line.split('\t')
        gsplit = gname.split('.')
        gene_id = '.'.join(gsplit[:len(gsplit)-1])
        annot_src = gsplit[-1]
        if annot_src not in annot_fhs:
            annot_fhs[annot_src] = open(annot_src+"."+annot_type+".sums.tsv","w")
            annot_fhs[annot_src].write("gene_id\ttotal_length\tchromosome\tstart\tend\t%s\n" % (counts_header))

        key = "\t".join([chrm, gstart, gend, strand])
        if key not in annot_map:
            annot_map[key] = []
        annot_map[key].append([gene_id, annot_src])

for line in sys.stdin:
    line = line.rstrip()
    fields = line.split('\t')
    (chrm, gstart, gend, strand)  = fields[:4]
    counts = "\t".join(fields[4:])
    #output is: Gene   bp_length       chromosome      start   end count1 count2 count3....
    outstr = "\t".join([str(int(gend) + int(gstart) + 1), chrm, gstart, gend])+"\t"+counts+"\n"
    key = "\t".join([chrm, gstart, gend, strand])
    for (gene_id, annot_src) in annot_map[key]:
        annot_fhs[annot_src].write(gene_id + "\t" + outstr)

for fh in annot_fhs.values():
    fh.close()
