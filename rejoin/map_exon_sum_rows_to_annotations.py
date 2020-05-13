import sys
from collections import Counter

#first load ',' delimited list of annotation suffixes
#e.g. "G026,G029,R109,ERCC,SIRV" 
annotations = sys.argv[2]
annotations = annotations.split(',')
num_annotations = len(annotations)
amap = {a:i for (i,a) in enumerate(annotations)}

#for masks, each annotation gets this number of digits
#this is so we can avoid having a delimiter in the masks string
num_digits_per_annotation = 3
if len(sys.argv) > 3:
    num_digits_per_annotation = int(sys.argv[3])

e2amap = {}
#second load exons2gene.annotation mapping
#e.g. G029.G026.R109.F006.20190220.gtf.exons2genes
with open(sys.argv[1],"r") as fin:
    for line in fin:
        #sys.stdout.write(line)
        fields = line.rstrip().split('\t')
        #use the exon_id|chr|start|end|strand to match
        #as exon IDs can be shared across annotations but not have the same exons (e.g. Gencode versions)
        eid = fields[0]
        annotation = None
        #special casing for the spike-ins
        first_4 = fields[1][:4]
        if first_4 == 'ERCC' or first_4 == 'SIRV':
            annotation = first_4
        #otherwise try to find the annotation in the string 
        #starting from the end
        #(where it normally is except for a minority of refseq genes)
        else:
            annot_suffixes = fields[1].split('.')
            i = len(annot_suffixes) - 1
            while i >= 0:
                if annot_suffixes[i] in amap:
                    break
                i -= 1
            if i < 0:
                continue
            annotation = annot_suffixes[i]
        if eid not in e2amap:
            e2amap[eid] = Counter()
        #update exon_id's count of its annotation's position in the count vector
        #e.g. for this exon final count vector could be:
        #"120300" meaning 1 row for G026, 2 copies of the exon row for G029, 
        #0 copies for R109, 3 copies for F006 and 0 copies each for ERCC, SIRV
        e2amap[eid][(amap[annotation])] += 1

#for each row, print N for each annotation which has the exon
#where N is >=1, determining the number of copies of the exon 
#row that should be printed in the final per-study-per-annotation file
for line in sys.stdin:
    fields = line.rstrip().split('\t')
    eids = fields[0].split(';')
    seen = Counter()
    mask = [0]*num_annotations
    #offset start by 1 so we're 1 base
    fields[2] = str(int(fields[2])+1)
    #chrm|start|end|strand
    e_suffix = '|'.join(['|'.join(fields[1:4]), fields[5]])
    for eid in eids:
        eid = eid + '|' + e_suffix
        if eid not in e2amap:
            continue
        seen.update(e2amap[eid])
        if len(seen.keys()) == num_annotations:
            break
    #now setup position count vector for this exon
    for apos in seen:
        mask[apos] = seen[apos]
    #print exactly num_annotations * num_digits_per_annotation number of digits for the mask
    mask_str = ''.join(['{0:0{1}.0f}'.format(m, num_digits_per_annotation) for m in mask])
    #dont skip writing out full 0 masks, we need to exactly match the existing exon sums file's
    #number of rows and they'll be filtered out later
    sys.stdout.write("%s\t%s" % (mask_str, line))
