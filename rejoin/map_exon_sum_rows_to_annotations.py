import sys

#first load ',' delimited list of annotation suffixes
#e.g. "G026,G029,R109,ERCC,SIRV" 
annotations = sys.argv[2]
annotations = annotations.split(',')
num_annotations = len(annotations)
amap = {a:i for (i,a) in enumerate(annotations)}

e2amap = {}
#second, load exonID2annotation mapping
with open(sys.argv[1],"r") as fin:
    for line in fin:
        #sys.stdout.write(line)
        fields = line.rstrip().split('\t')
        eid = fields[0].split('|')[0]
        annotation = None
        first_4 = fields[1][:4]
        #special casing for the spike-ins
        if first_4 == 'ERCC' or first_4 == 'SIRV':
            annotation = first_4
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
            e2amap[eid] = set()
        e2amap[eid].add(amap[annotation])

#for each row, add the int IDs of the annotation as a bitmask
for line in sys.stdin:
    fields = line.rstrip().split('\t')
    eids = fields[0].split(';')
    seen = set()
    mask = ["0"]*num_annotations
    for eid in eids:
        if eid not in e2amap:
            continue
        seen = seen.union(e2amap[eid])
        if len(seen) == num_annotations:
            break
    for abit in seen:
        mask[abit] = "1"
    mask_str = ''.join(mask)
    if "1" not in mask_str:
        continue
    sys.stdout.write("%s\t%s" % (mask_str, line))
