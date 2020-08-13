#!/usr/bin/env python3.6
import sys

QC_START_COL=70
QC_END_COL=129

new_header_fields = ["star.number_of_input_reads_both",
"star.all_mapped_reads_both",
"star.number_of_chimeric_reads_both",
"star.number_of_reads_mapped_to_multiple_loci_both",
"star.number_of_reads_mapped_to_too_many_loci_both",
"star.number_of_reads_unmapped:_other_both",
"star.number_of_reads_unmapped:_too_many_mismatches_both",
"star.number_of_reads_unmapped:_too_short_both",
"star.uniquely_mapped_reads_number_both",
"star.%_mapped_reads_both",
"star.%_chimeric_reads_both",
"star.%_reads_mapped_to_multiple_loci_both",
"star.%_reads_mapped_to_too_many_loci_both",
"star.%_reads_unmapped:_other_both",
"star.%_reads_unmapped:_too_many_mismatches_both",
"star.%_reads_unmapped:_too_short_both",
"star.uniquely_mapped_reads_%_both"]

total_col = 0
fields_to_add=[102,82,100,104,106,108,110,112,128]

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
        for i in fields_to_add:
            f1 = float(fields[i])
            f2 = float(fields[i+1])
            new_qc[j] = f1 + f2
            if(j > 0):
                new_qc_pct[j-1] = '%.1f' % (100*(new_qc[j] / new_qc[total_col]))
            j += 1
        new_qc_str = '\t'.join([str(x) for x in new_qc]) + '\t' + '\t'.join([str(x) for x in new_qc_pct])
        sys.stdout.write(line + '\t' + new_qc_str + '\n')
