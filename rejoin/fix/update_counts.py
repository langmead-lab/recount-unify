#!/usr/bin/env python3
import sys

#e.g. gene.counts.G026
#first column is gene_id, rest are sums
#NO HEADER
updated_counts_for_annotationF=sys.argv[1]
sanity_check_geneF=sys.argv[2]
#sys.stderr.write(sanity_check_gene+"\n")
#only expect 1 line
with open(sanity_check_geneF,"r") as fin0:
    for line in fin0:
        (chrm,start,end,strand,sanity_gene)=line.rstrip().split('|')

new_counts_map = {}
with open(updated_counts_for_annotationF,"r") as fin:
    for line in fin:
        line = line.rstrip()
        gene_id_pos = line.find("\t")
        if gene_id_pos == -1:
            sys.stderr.write("bad row, no tabs, terminating early!")
            sys.exit(-1)
        gene_id = line[:gene_id_pos]
        new_counts_map[gene_id] = line[gene_id_pos:]

#has 3 header lines
idx=0
for line in sys.stdin:
    line_original = line
    idx += 1
    if idx < 4:
        sys.stdout.write(line)
        continue
    line = line.rstrip()
    gene_id_pos = line.find("\t")
    if gene_id_pos == -1:
        sys.stderr.write("bad row, no tabs, terminating early!")
        sys.exit(-1)
    current_counts = line[gene_id_pos:]
    gene_id = line[:gene_id_pos]
    if gene_id not in new_counts_map:
        sys.stdout.write(line_original)
        continue
    new_counts = new_counts_map[gene_id]
    if new_counts == current_counts:
        sys.stderr.write("gene_counts_match\t"+gene_id+"\n")
        sys.stdout.write(line_original)
    else:
        if gene_id == sanity_gene:
            sys.stderr.write("sanity_check_gene_counts_changed\t"+sanity_gene+"\n")
            sys.exit(-1)
        sys.stderr.write("gene_counts_no_match\t"+gene_id+"\n")
        sys.stdout.write(gene_id+new_counts+"\n")
