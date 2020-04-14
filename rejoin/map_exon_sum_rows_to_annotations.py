import sys

#',' delimited list of annotation suffixes
#e.g. "G026,G029,R109,ERCC,SIRV" 
annotations = sys.argv[2]
annotations = annotations.split(',')
num_annotations = len(annotations)
amap = {a:i for (i,a) in enumerate(annotations)}

e2amap = {}
#first, load exonID2annotation mapping
with open(sys.argv[1],"r") as fin:
    for line in fin:
        #sys.stdout.write(line)
        fields = line.rstrip().split('\t')
        eid = fields[0]
        annotations = fields[1].split(',')
        annot_len = 0
        for annot in annotations:
            if eid not in e2amap:
                e2amap[eid] = set()
            e2amap[eid].add(amap[annot])

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
    sys.stdout.write("%s\t%s" % (''.join(mask), line))
