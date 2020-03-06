import sys
sums=[]
sum_len=0
count=0
#expected timings columns in input
#total_runtime   align   align_unmapped  all     bamcount        bamcount_unmapped       download        exon_fc_count_all       exon_fc_count_unique      extract_jx      fastq_check     gene_fc_count_all       gene_fc_count_unique    make_manifest   remove_fastqs   remove_unmapped_fastqs    salmon  sort
min_time = 0
max_time = 0
for line in sys.stdin:
    count += 1
    if sum_len == 0:
        sums = [float(z) for z in line.rstrip().split('\t')]
        sum_len = len(sums)
        min_time = sum(sums[1:])
        max_time = min_time
        continue
    sums_ = [float(z) for z in line.rstrip().split('\t')]
    min_time = min([sum(sums_[1:]), min_time])
    max_time = max([sum(sums_[1:]), max_time])
    sums = [sums[i]+z for (i,z) in enumerate(sums_)]
parallel_total_runtime = sum(sums[1:])
sys.stdout.write('totals:\t%s\t%d\t%s\n' % (str(parallel_total_runtime), count,'\t'.join([str(z) for z in sums])))
align = round(100*sums[1]/parallel_total_runtime,2)
download = round(100*sums[6]/parallel_total_runtime,2)
bamcount = round(100*sums[4]/parallel_total_runtime,2)
every_but_align_download = 100-(download+align)
time_per_run = parallel_total_runtime/count
sys.stdout.write("percents:\talignment\tdownload\teverything_else\tavg_time_per_run\tmin_time_per_run\tmax_time_per_run\n")
sys.stdout.write("percents:\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n" % (align, download, every_but_align_download, round(time_per_run), round(min_time), round(max_time)))
sys.stdout.write("bamcount_pct:%s\n" % bamcount)
