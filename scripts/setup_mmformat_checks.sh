d=$1
pushd $d/study
ls *.MM.gz > all_mms
cat  all_mms | perl -ne 'chomp; $f=$_; $f2=$f; $f2=~s/\.MM\.gz/\.RR\.gz/; print "../../scripts/check_mmformat_output.sh $f $f2 sra\n";' > all_mms.check.jobs
parallel -j 30 < all_mms.check.jobs > all_mms.check.jobs.run 2>&1
fgrep "gz	OK" all_mms.check.jobs.run | fgrep -v "echo" | cut -f 1 > all_mms.check.jobs.run.good
fgrep -v -f all_mms.check.jobs.run.good all_mms > all_mms.bad
for f in `cat all_mms.bad`; do
    mv $f ${f}.empty
done
ls -tlr *.empty
popd
