import sys
line=sys.argv[1]
complement=len(sys.argv) > 2
#find contiguous ranges of 
#sample positions to condense
#this keeps the command line smaller
prev_run_pos = -100
start_pos = -100
positions = []

if complement:
    bad = {}
    start_pos = 1
    for run_pos in line.split(','):
        end_pos = run_pos - 1
        #skipping position, can at least happen at position 1
        if end_pos < start_pos:
            start_pos = run_pos + 1
            continue
        #only including one position
        elif start_pos == end_pos:
            positions.append(str(start_pos))
            start_pos = run_pos + 1
        else:
            pos = "%d-%d" % (start_pos, end_pos)
            positions.append(pos)
            start_pos = run_pos + 1
    pos = "%d-" % (start_pos)
    positions.append(pos)
    cut_positions = ','.join(positions) 
    sys.stdout.write("cut -f"+cut_positions)
    sys.exit(0)

#positions are what we want to keep
for run_pos in line.split(','):
    if run_pos != prev_run_pos+1:
        pos = prev_run_pos
        if start_pos != pos:
            pos = "%d-%d" % (start_pos, prev_run_pos)
        if prev_run_pos != -100:
            positions.append(str(pos))
        start_pos = run_pos 
    prev_run_pos = run_pos
pos = prev_run_pos
if start_pos != pos:
    pos = "%d-%d" % (start_pos, prev_run_pos)
if prev_run_pos != -100:
    positions.append(str(pos))
cut_positions = ','.join(positions) 
sys.stdout.write("cut -f"+cut_positions)
