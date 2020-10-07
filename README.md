# recount-unify
Next step after recount-pump in the Monorail pipeline

This step summarizes the gene, exon, and junction level counts which are produced by the `recount-pump` workflow on a per-sequencing run basis.

Gene and exon level counts are produced as full matrices per annotation, per study (if running on an SRA tranche of multiple studies and with multiple annotations, e.g. GencodeV26 & RefSeq). 

Junctions are produced currently as per the Snaptron format, i.e. multiple studies are kept in the same Snaptron compilation if doing SRA, also the annotations used to mark junctions as annotated (or novel) is a separate, and typically larger set, than the annotations used in the gene/exon counts.

## Dependencies

* `Snakemake` is needed for running unifying workflow
* `python3` is used to run the main Snakemake workflow (though `python2` will probably work)
* `pypy` is needed for certain steps within the workflow (jx merging and annotating)
* `zstd` is needed for decompressing `recount-pump` outputs
* `pigz` is needed for compressing final outputs

While the Snakemake file does most of the heavy lifting, there are a number of separate steps which are still outside the Snakemake workflow.

The following files/information are needed to run the unifier, but are specific to your project:

* a tab-delimited file of the project study(s) and sample/run ids used in the `recount-pump` run (e.g. `project1.tsv` below)
* compilation ID for your project that doesn't collide with existing recount3/Snaptron2 compilations IDs (`compilation_id`)
* number of samples/runs in the project (`#_samples` below, should be exactly the set of runs/samples which successfully ran through `recount-pump`)

The following are more generic files used across projects sharing the same reference genome (e.g. `hg38`) and set of annotations, e.g. `G026,G029,R109,ERCC,SIRV,F006`:

* set of annotation short names used (e.g. `G026,G029,R109,ERCC,SIRV,F006`)
* disjoint exon-to-gene mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed`)
* disjoint exon-to-annotated exon mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed`)
* final gene disjoint-to-annotation mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed`)
* a list of annotated junctions (e.g. `annotated_junctions.tsv.gz`)
* the list of annotated exon intervals used in the `recount-pump` run (e.g. `exons.bed.w_header.gz`)
  this file has to be be exactly the same order as the exon sums that's produced by `bamcount` as part of `recount-pump`
* genome reference chromosome sizes file (e.g. `hg38.recount_pump.fa.sizes.tsv`)

These files can all be downloaded from https://github.com/langmead-lab/monorail-run-data
for recount3 human & mouse projects.

And finally, also needed is the genome reference chromosome FASTA file (e.g. `hg38.recount_pump.fa`), which needs to be exactly the same as what's used in `recount-pump`, available from https://recount-ref.s3.amazonaws.com/hg38/fasta.tar.gz
The chromosome names need to exactly match what's in genome reference chromosome sizes file above, or this file should be regenerated based on the FASTA file.

## Prep

This step assumes that all checking/filtering of `recount-pump` output has already been done previously.

First, you need to setup the unifier working subdirectory (of the unifier github checkout):

```mkdir project1
cd project1
python ../sample_ids/assign_compilation_ids.py --accessions-file project1.tsv --compilation-code <compilation_code> > project1.ids.tsv
```

Where `<compilation_code>` needs to checked against the list of existing recount3/Snaptron2 compilations so it doesn't collide, https://github.com/langmead-lab/monorail-run-data/blob/master/recount-unify/sample_ids/all_compilation_codes.tsv

Then, you need to create symlinks from the output of `recount-pump` (assumes the same filesystem):

```/bin/bash -x ../scripts/find_done.sh /path/to/project_recount-pump_output links project_recount-pump_file_prefix```

An example for an human SRA tranche (5) on the IDIES systems:

```/bin/bash -x ../scripts/find_done.sh /full/path/to/langmead/sra_human_v3_5 links sra_human_v3```

Finally, a file named `blank_exon_sums` needs to be in the working directory of the unifier.

This file should contain the exact number of lines that are in the list of annotated exon intervals used in the `recount-pump` run sans the header line (e.g. `exons.bed.w_header.gz`) above, but where every line is a single 0.  This is used as a placeholder for samples which for some reason don't have sums, as the total # of columns in sums output needs to match the total # of input samples.

The blank_exon_sums files in the root of this repo was used for SRAv3/GTExv2/TCGAv2 (human).

## Layout of links to recount-pump output

Due to the importance of this part, this get its own section.

The `scripts/find_done.sh` script mentioned in the previous section *should* organize the symlinks to the original recount-pump output directories correctly, however, it's worth checking given that the rest of the Unifier is critically sensitive to how the links are organized.

For example, if you find that you're getting blanks instead of actual integers in the `all.exon_bw_count.pasted.gz` file, it's likely a sign that the input directory hierarchy was not laid out correctly.

Assuming your top level directory for input is called `links`, the expected directory hierarchy for each sequencing run/sample is:

`links/<study_loworder>/<study>/<run_loworder>/<run>/<symlink_to_recount-pump_attempt_director_for_this_run>

e.g.:

`links/94/SRP019994/83/SRR797083/sra_human_v3_41_in26354_att2`

where `sra_human_v3_41_in26354_att2` is the symlink to the actual recount-pump generated attempt director for run `SRR797083` in study `SRP019994`.

`study_loworder` and `run_loworder` are *always* the last 2 characters of the study and run accessions/IDs respectively.

## Running the Unifier Workflow

An semi-generic example of the unifier command for a human SRA tranche:

```
snakemake -j <#_threads> --stats ./stats.json --snakefile ../Snakefile -p --config input=links staging=unified sample_ids_file=project1.ids.tsv annotated_sjs=/path/to/annotated_junctions.tsv.gz existing_sums=/path/to/exons.bed.w_header.gz compilation_id=<compilation_id> gene_rejoin_mapping=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed exon_rejoin_mapping=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed num_samples=<#_samples> ref_sizes=/path/to/hg38.recount_pump.fa.new_sizes ref_fasta=/path/to/hg38.recount_pump.fa recount_pump_output=/path/to/project_recount-pump_output gene_mapping_final=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed annotation_list=G026,G029,R109,ERCC,SIRV,F006
```

`links` is the directory created in the prep step above which contains the symlinks to to all of recount-pump's ourput files.

`unified` is a directory which will be created by the Snakemake workflow which contains intermediate files, once a run's final files have been produced and checked this directory can be removed.


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
