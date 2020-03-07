import sys
from collections import Counter

samples2cov = Counter()
for line in sys.stdin:
	fields = line.rstrip().split('\t')
	for combo in fields[11].split(',')[1:]:
		(sid,cov)=combo.split(':')
		samples2cov[sid] += int(cov)

for sample in samples2cov.keys():
	sys.stdout.write("%s\t%d\n" % (sample, samples2cov[sample]))
