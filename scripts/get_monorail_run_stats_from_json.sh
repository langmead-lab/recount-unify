
#to be run in the post-recount-pump output for the project/tranche
#should include 1) spike-ins 2) overlapping runs from other tranches for completeness

dir=$(dirname $0)

#the manifest file is the signal that the job completed (if it's > 1000 bytes)
find . -name "*.manifest" -size +1000c -exec ls -l {} \; > all_non0_manifests
cat all_non0_manifests | tr -s " " \\t | cut -f 9 | sed 's/sra.manifest/sra.stats.json/' > all_stats_jsons
cat all_stats_jsons | perl -ne 'chomp; print "pypy '$dir'/parse_snakemake_stats.py $_ > $_.parsed 2> $_.parsed.err\n";' > stats_json.parse.jobs
