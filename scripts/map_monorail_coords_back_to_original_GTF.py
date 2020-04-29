import sys
import re

coord_map = {}
with open(sys.argv[1],"r") as fin:
    for (idx,line) in enumerate(fin):
        (gid, bp_length, chrm, start, end) = line.rstrip().split('\t')
        #coord_map['_'.join([gid,chrm,start,end])] = [bp_length, chrm, start, end]
        coord_map['_'.join([gid,chrm,start,end])] = [idx,bp_length]

#exon or gene
feat_type = sys.argv[2]
#4 char. abbrev for annotation, e.g. R109 for refseq hg38
annot = sys.argv[3]

#GL000008.2^ICurated Genomic^Iexon^I124376^I125329^I.^I-^I.^Itranscript_id "gene14440.R109"; gene_id "gene14440.R109"; gene_name "SNX18P15"; Dbxref "GeneID:100419019,HGNC:HGNC:39623"; Name "SNX18P15"; description "sorting nexin 18 pseudogene 15"; gbkey "Gene"; gene_biotype "pseudogene"; pseudo "true";
patt = re.compile(r'gene_id\s+"([^"]+)"')
for line in sys.stdin:
    if line[0] == '#':
        sys.stdout.write(line) 
        continue
    fields = line.rstrip().split('\t')
    if fields[2] != feat_type:
        continue
    (chrm, start, end) = (fields[0], fields[3], fields[4])
    fields[8] = fields[8].replace('.%s' % (annot),'')
    m = patt.search(fields[8])
    if m is None:
        sys.stderr.write("MISSING_gene_id, skipping %s" % (line))
        continue
    gid = m.group(1)
    #gid = gid.replace('.%s' % (annot),'',1)
    key = '_'.join([gid,chrm,start,end])
    if key not in coord_map:
        sys.stderr.write("NOT_USED_IN_MONORAIL_gene_id, skipping %s" % (line))
        continue
    fields[5] = coord_map[key][1]
    #output a key we can sort later on to match the exact order of the sums file
    sys.stdout.write(str(coord_map[key][0])+'\t'+'\t'.join(fields)+'\n') 
