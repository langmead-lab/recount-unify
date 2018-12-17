# recount-unify
next step after recount-pump in the monorail pipeline

## To run the patroller script:
`python3 ./scripts/patroller.py --num-sm-proc 8 --snakefile ./Snakefile --annotated-sj-file /path/to/annotated_junctions.tsv.gz --sample-ID-file /path/to/samples.tsv --incoming-dir /path/to/recount_pump_output --existing-exon-sums-file exons.bed.gz --debug --project <project_name>`

where `--num-sm-proc` is the number of concurrent Snakemake processes you want to run (i.e. `snakemake -j`),
`--existing-exon-sums-file` is a tab delimited list of just chromosome,start,ends,
`<project_name>` is for downloading the project study2run mapping from S3 (e.g. `srav1`).


merge is for merging junctions from all samples into one sparse matrix

disjoin is for taking an annotation and creating a disjoint set of exons from it
it also contains code to re-assemble counts for transcripts, genes, original exons from the disjoint exon counts

annotate/annotation contains scripts for compiling a multi-source junction annotation

annotate contains script(s) for annotating the merged junction file
