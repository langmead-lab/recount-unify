#!/usr/bin/env python3
#ensure that the order of gene counts for a particular annotation
#matches what's in the recount3 GTF/GFF file, otherwise the counts won't load
#could use a lot of RAM as it has to store the gene counts matrix in memory
import sys

#use pigz to do the decompression and redirect to this argument
#e.g. <(pcat $gene_countsF)
gene_countsF=sys.argv[1]

counts_map={}
#has 3 header lines
idx=0
with open(gene_countsF,"r") as fin:
    for line in fin:
        idx += 1
        if idx < 4:
            sys.stdout.write(line)
            continue
        gene_id_pos = line.find("\t")
        if gene_id_pos == -1:
            sys.stderr.write("bad row, no tabs, terminating early!")
            sys.exit(-1)
        gene_id = line[:gene_id_pos]
        counts_map[gene_id] = line[gene_id_pos:]

for line in sys.stdin:
    line = line.rstrip()
    counts = counts_map[line]
    sys.stdout.write(line+counts)
    del counts_map[line]

if len(counts_map) > 0:
    sys.stderr.write("left over gene counts, terminating early!\n")
    sys.exit(-1)
