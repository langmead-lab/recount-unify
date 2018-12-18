# recount-unify
next step after recount-pump in the monorail pipeline

## To run the patroller script:
`python3 ./patroller.py --num-sm-proc #_procs --snakefile ./Snakefile --annotated-sj-file /path/to/annotated_junctions.tsv.gz --sample-ID-file /path/to/samples.tsv --incoming-dir /path/to/recount_pump_output --staging-dir /path/to/dir_to_hardlink_into --output-dir /path/to/dir_to_store_final_unified_output --existing-exon-sums-file /path/to/exons.bed.gz --project <project_name> --debug`

where:
* `--num-sm-proc` is the number of concurrent Snakemake processes you want to run (i.e. `snakemake -j`)
* `--annotated-sj-file` takes a gzipped, tab-delimited list of known annotated junctions and their source annotation abbreviations
* `--sample-ID-file` takes a uncompressed, tab-delimited file of rail_id2run_accession mappings
* `--incoming-dir` takes the top-level directory where all recount-pump output is dumped
* `--staging-dir` takes a path to store all the hardlinks to all files in the recount_pump_output dir
* `--output-dir ` takes a path to store the output of the unifier process which will be organized along the same lines as the --incoming-dir
* `--existing-exon-sums-file` is a gzipped, tab-delimited list of just chromosome,start,ends
*`<project_name>` is for downloading the project study2run mapping from S3 (e.g. `srav1`)

`--debug` will force the pipeline to only run upto 5 random studies and then quit, otherwise it'll do every study that
has a complete set of "done" runs and then loop forever after sleeping for a few seconds at each loop.

### Patroller Dependencies
* zstd executable in the path
* Snakemake executable in the path
* Python 3.6
* Python 2.7

## Other functionality

*merge* is for merging junctions from all samples into one sparse matrix

*disjoin* is for taking an annotation and creating a disjoint set of exons from it
it also contains code to re-assemble counts for transcripts, genes, original exons from the disjoint exon counts

*annotate/annotation* contains scripts for compiling a multi-source junction annotation

annotate contains script(s) for annotating the merged junction file
