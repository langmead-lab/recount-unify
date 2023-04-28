#!/usr/bin/env python3.6
import sys

new_header_fields = ["star.number_of_input_reads_both",
"star.all_mapped_reads_both",
"star.number_of_chimeric_reads_both",
"star.number_of_reads_mapped_to_multiple_loci_both",
"star.number_of_reads_mapped_to_too_many_loci_both",
"star.uniquely_mapped_reads_number_both",
"star.%_mapped_reads_both",
"star.%_chimeric_reads_both",
"star.%_reads_mapped_to_multiple_loci_both",
"star.%_reads_mapped_to_too_many_loci_both",
"star.uniquely_mapped_reads_%_both"]

total_col = 0
add_both = True

#just take original QC stats output file which just have run,study as extra columns
#fields_to_add=[63,43,61,65,67,69,71,73,89]
#fields_to_add=[64,44,62,66,68,90]
#fields_to_add=[50,31,48,52,54,70]
fields_to_add=[51,32,49,53,55,71]
add_both = False

with open(sys.argv[1],"r") as fin:
    header_fields = []
    hlen = 0
    flen = len(fields_to_add)
    new_qc = [0.0]*flen
    new_qc_pct = ['0.0']*(flen-1)
    for line in fin:
        line = line.rstrip()
        fields = line.split('\t')
        if hlen == 0:
            header_fields = fields
            hlen = len(fields)
            sys.stdout.write(line + '\t' + '\t'.join(new_header_fields) + '\n')
            continue
        j = 0
        bad_row = False
        for i in fields_to_add:
            f1 = float(fields[i])
            new_qc[j] = f1
            if add_both:
                f2 = float(fields[i+1])
                new_qc[j] = f1 + f2
            if j == 0 and new_qc[j] == 0.0:
                new_qc[total_col]=0.000000001
                bad_row = True
            if(j > 0):
                new_qc_pct[j-1] = '%.1f' % (100*(new_qc[j] / new_qc[total_col]))
            j += 1
        if bad_row:
            sys.stderr.write("BAD_ROW_0_TOTAL_INPUT_READS\t"+line+'\n')
        new_qc_str = '\t'.join([str(x) for x in new_qc]) + '\t' + '\t'.join([str(x) for x in new_qc_pct])
        sys.stdout.write(line + '\t' + new_qc_str + '\n')
