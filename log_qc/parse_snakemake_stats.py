import json
import sys

#default runtime to get
dt='min-runtime'
fname = sys.argv[1]
basename = fname.split('/')[-1]
fnames = basename.split('!')
study=fnames[1]
run=fnames[0]
fh = open(fname,"rb")
ds = json.load(fh)
output = []
total = ds['total_runtime']
ds = ds['rules']
#sys.stderr.write("study\trun\ttotal_runtime\t%s\n" % ('\t'.join([rule for rule in sorted(ds.keys())])))
sys.stdout.write("%s\t%s\t%s\t%s\n" % (study,run,str(total),'\t'.join([str(ds[rule][dt]) for rule in sorted(ds.keys())])))	
#download = ds['download'][dt]
#alignment = ds['alignment'][dt]
