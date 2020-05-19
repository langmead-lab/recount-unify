# recount-unify
Next step after recount-pump in the monorail pipeline

This step summarizes the gene, exon, and junction level counts which are produced by the Monorail pipeline on a per-sequencing run basis.

Gene and exon level counts are produced as full matrices per annotation, per study (if running on an SRA tranche of multiple studies and with multiple annotations, e.g. GencodeV26 & RefSeq). 

Junctions are produced currently as per the Snaptron format, i.e. multiple studies are kept in the same Snaptron compilation if doing SRA, also the annotations used to mark junctions as annotated (or novel) is a separate, and typically larger set, than the annotations used in the gene/exon counts.

## Dependencies

* Snakemake is needed for running unifying workflow
* python3 is used to run the main Snakemake workflow (though python2 will probably work)
* pypy is needed for certain steps within the workflow (jx merging and annotating)
* zstd is needed for decompressing Monorail outputs
* pigz is needed for compressing final outputs

While the Snakemake file does most of the heavy lifting, there are a number of separate steps which are still outside the Snakemake workflow.

The following files/information are needed to run the unifier, but our specific to your project:

* a tab-delimited file of the project study(s) and sample/run ids used in the `recount-pump` run (e.g. `project1.ids.tsv` below)
* compilation ID for your project that doesn't collide with existing recount3/Snaptron2 compilations IDs (`compilation_id`)
* number of samples/runs in the project (`#_samples` below, should be exactly the set of runs/samples which successfully ran through `recount-pump`)

The following are more generic files used across projects sharing the same reference genome (e.g. `hg38`) and set of annotations, e.g. `G026,G029,R109,ERCC,SIRV,F006`:

* a list of annotated junctions (e.g. `annotated_junctions.tsv.gz`)
* disjoint exon-to-gene mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed`)
* disjoint exon-to-annotated exon file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed`)
* final gene disjoint-to-annotation mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed`)
* set of annotation short names used (e.g. `G026,G029,R109,ERCC,SIRV,F006`)
* genome reference chromosome sizes file (e.g. `hg38.recount_pump.fa.new_sizes`)
* genome reference chromosome FASTA file (e.g. `hg38.recount_pump.fa`), needs to be exactly the same as what's used in `recount-pump`
* the list of annotated exon intervals used in the `recount-pump` run (e.g. `exons.bed.w_header.gz`)
  this file has to be be exactly the same order as the exon sums that's produced by `bamcount` as part of `recount-pump`

These files can all be downloaded from https://github.com/langmead-lab/monorail-run-data
for recount3 human & mouse projects.

## Prep

This step assumes that all checking/filtering of Monorail output has already been done previously.

First, you need to setup the unifier working subdirectory (of the unifier github checkout):

```mkdir project1
cd project1
python ../sample_ids/assign_compilation_ids.py --accessions-file project1.tsv --compilation-code <compilation_code> > project1.ids.tsv
```

Where `<compilation_code>` needs to checked against the list of existing recount3/Snaptron2 compilations so it doesn't collide.

Then, you need to create symlinks from the output of Monorail (assumes the same filesystem):

```/bin/bash -x ../scripts/find_done.sh /path/to/project_monorail_output links project_monorail_file_prefix```

An example for an human SRA tranche (5) on the IDIES systems:

```/bin/bash -x ../scripts/find_done.sh /scratch/ceph/langmead/sra_human_v3_5 links sra_human_v3```

## Running the Unifier Workflow

An example of the unifier command for human SRA tranche 1 on IDIES:

```
snakemake -j <#_threads> --stats ./stats.json --snakefile ../Snakefile -p --config input=links staging=unified sample_ids_file=project1.ids.tsv annotated_sjs=/path/to/annotated_junctions.tsv.gz existing_sums=/path/to/exons.bed.w_header.gz compilation_id=<compilation_id> gene_rejoin_mapping=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed exon_rejoin_mapping=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed num_samples=<#_samples> ref_sizes=/path/to/hg38.recount_pump.fa.new_sizes ref_fasta=/path/to/hg38.recount_pump.fa recount_pump_output=/path/to/project_monorail_output gene_mapping_final=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed annotation_list=G026,G029,R109,ERCC,SIRV,F006
```


## QC summary

This is typically run after the unifier has been run on the project/tranche (if you want it to include the intron sums) in the unifier working subdirectory of the project/tranche:

```python3 ../log_qc/parse_logs_for_qc.py --incoming-dir links --sample-mapping ids.tsv --intron-sums intron_counts_summed.tsv > qc.tsv 2> qc.err```


## [DEPRECATED] To run the patroller script:
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

Before running, make sure you've created the `--staging-dir`, the `--output-dir`, and a subdirectory in the working directory `logs` where per-study output from the Snakemake will be deposited.  Make the staging, output and logs directories be subdirectories of `$PWD` when running script.  Assumption is that `$PWD/logs` subdirectory exists.

### Patroller Dependencies
* zstd executable in the path
* Snakemake executable in the path
* Python 3.6
* pypy 2.7

## Other functionality

*merge* is for merging junctions from all samples into one sparse matrix

*disjoin* is for taking an annotation and creating a disjoint set of exons from it
it also contains code to re-assemble counts for transcripts, genes, original exons from the disjoint exon counts

*annotate/annotation* contains scripts for compiling a multi-source junction annotation

annotate contains script(s) for annotating the merged junction file
