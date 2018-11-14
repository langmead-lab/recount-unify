python2 merge.py --list-file test/test_samples.tsv  > test_samples.tsv.results
diff test/true1abc.tsv test_samples.tsv.results

python2 merge.py --list-file test/test_samples2.tsv --append-samples --existing-jx-db test/test2.existing > test_samples2.tsv.results
diff test/true2_9.tsv test_samples2.tsv.results
