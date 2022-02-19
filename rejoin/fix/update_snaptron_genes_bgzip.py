#!/usr/bin/env python3
import sys
#run like:
#/usr/bin/time -v pcat genes.bgz.pre_rejoin_fix 2> genes.updated.run | python3 ../update_snaptron_genes_bgzip.py gtex.studies.snaptron.pasted.recalulated_stats2 > genes.updated &

#then, if successful, run this to update the sqlite3 db:
#cat gene_updatesF.wsnapids | perl -ne 'chomp; $f=$_; ($snap_id,$gid,$samples,$slen,$ssum,$savg,$smedian,$src)=split(/\t/,$f,-1); print "UPDATE intron SET samples=\"$samples\",samples_count=$slen,coverage_sum=$ssum,coverage_avg=$savg,coverage_median=$smedian,source_dataset_id=$src WHERE snaptron_id=$snap_id;\n";' > genes.updates.sql
#takes ~21s on snaptron01 on nvme2 to do gtexv2:
#/usr/bin/time -v sqlite3 genes.sqlite3 < genes.updates.sql > genes.updates.sql.run 2>&1

gene_updatesF=sys.argv[1]
update_map = {}
with open(gene_updatesF,"r") as fin:
    for line in fin:
        fields = line.split('\t')
        gid = fields[0]
        newline = '\t'.join(fields[1:])
        update_map[gid] = newline

wsnap_ids = open("gene_updatesF"+".wsnapids","w")
for line in sys.stdin:
    first_comma_pos = line.find(',')
    #if there's no sample list, we just write it out
    if first_comma_pos == -1:
        sys.stdout.write(line)
        continue
    #split up to, but not incuding tab before first comma
    ffields = line[:first_comma_pos-1].split('\t')
    if len(ffields) != 11:
        sys.stderr.write("bad comma position seek for " + line)
        sys.exit(-1)
    prefix = line[:first_comma_pos-1]
    gene_tab_pos = prefix.rfind('\t')
    gene_colon_pos = prefix.find(':')
    gene_id = prefix[gene_tab_pos+1:gene_colon_pos]

    if gene_id not in update_map:
        sys.stdout.write(line)
        continue
    newline = update_map[gene_id]
    wsnap_ids.write(ffields[0] + '\t' + gene_id + '\t' + newline)
    #replace previous samples + summary stats + datasource id with new version 
    sys.stdout.write(prefix + '\t' + newline) 

wsnap_ids.close()
