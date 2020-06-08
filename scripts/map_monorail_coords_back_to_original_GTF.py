import sys
import re

#exon or gene
feat_type = sys.argv[2]
#4 char. abbrev for annotation, e.g. R109 for refseq hg38
annot = sys.argv[3]

delim = '\t'
if feat_type == 'exon':
    delim = '|'

coord_map = {}
with open(sys.argv[1],"r") as fin:
    for (idx,line) in enumerate(fin):
        if line[0] == '#' or 'chromosome' in line:
            continue
        fields = line.rstrip().split(delim)
        if feat_type == 'gene':
            (gid, bp_length, chrm, start, end) = fields
            coord_map[delim.join([gid,chrm,start,end])] = [idx,bp_length]
            continue
        (chrm, start, end, strand) = fields
        bp_length = str((int(end) - int(start)) + 1)
        key = delim.join([chrm,start,end,strand])
        #sys.stderr.write(key+'\n')
        coord_map[key] = [idx,bp_length]


#GL000008.2^ICurated Genomic^Iexon^I124376^I125329^I.^I-^I.^Itranscript_id "gene14440.R109"; gene_id "gene14440.R109"; gene_name "SNX18P15"; Dbxref "GeneID:100419019,HGNC:HGNC:39623"; Name "SNX18P15"; description "sorting nexin 18 pseudogene 15"; gbkey "Gene"; gene_biotype "pseudogene"; pseudo "true";
gpatt = re.compile(r'gene_id\s+"([^"]+)"')
for line in sys.stdin:
    if line[0] == '#':
        sys.stdout.write(line) 
        continue
    fields = line.rstrip().split('\t')
    if fields[2] != feat_type:
        continue
    (chrm, start, end, strand) = (fields[0], fields[3], fields[4], fields[6])
    fields[8] = fields[8].replace('.%s' % (annot),'')
    key = delim.join([chrm,start,end,strand])
    if feat_type == 'gene':
        m = gpatt.search(fields[8])
        if m is None:
            sys.stderr.write("MISSING_%s_id, skipping %s" % (feat_type, line))
            continue
        gid = m.group(1)
        key = delim.join([gid,chrm,start,end])
        #some refseq genes have their chromosome as a suffix
        if key not in coord_map:
            key = delim.join([gid+'.%s'%(chrm),chrm,start,end])
    if key not in coord_map:
        sys.stderr.write("NOT_USED_IN_MONORAIL_%s, key=%s, skipping %s" % (feat_type, key, line))
        continue
    fields[5] = coord_map[key][1]
    recount_exon_id = ''
    if feat_type == 'exon':
        recount_exon_id = ' recount_exon_id "%s";' % (key)
    #output a key we can sort later on to match the exact order of the sums file
    sys.stdout.write(str(coord_map[key][0])+'\t'+'\t'.join(fields)+recount_exon_id+'\n') 
