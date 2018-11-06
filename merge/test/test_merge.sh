python2 merge.py --list-file test/test_samples.tsv  > test_samples.tsv.results
diff test/true1abc.tsv test_samples.tsv.results
