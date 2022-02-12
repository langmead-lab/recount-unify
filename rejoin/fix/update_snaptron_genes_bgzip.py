#!/usr/bin/env python3
import sys

gene_updatesF=sys.argv[1]
update_map = {}
with open(gene_updatesF,"r") as fin:
    for line in fin:
        fields = line.split('\t')
        gid = fields[0]
        newline = '\t'.join(fields[1:])
        update_map[gid] = newline

for line in sys.stdin:
    first_comma_pos = line.find(',')
    #split up to, but not incuding tab before first comma
    ffields = line[:first_comma_pos-1].split('\t')
    if len(ffields) != 11:
        sys.stderr.write("bad comma position seek for " + line)
        sys.exit(-1)
    prefix = line[:first_comma_pos-1]
    gene_tab_pos = prefix.rfind('\t')
    gene_colon_pos = prefix.find(':')
    gene_id = prefix[gene_tab_pos+1:gene_colon_pos]

    newline = update_map[gene_id]
    #replace previous samples + summary stats + datasource id with new version 
    sys.stderr.write(prefix + '\t' + newline) 
