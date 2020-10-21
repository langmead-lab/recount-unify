# recount-unify
These instructions are primarily aimed at those who are internal to the groups involved in maintaining recount and Snaptron, not for most external users.
If you are an external user (not in that group), please start with:

https://github.com/langmead-lab/monorail-external#getting-reference-indexes

and

https://github.com/langmead-lab/monorail-external#unifier-aggregation-over-per-sample-pump-outputs


This phase of Monorail has two parts (2 Snakemake files) which both aggregate across the per-sample outputs produced by the eariler, `recount-pump` workflow phase of Monorail:

* base coverage summed for exons and genes based on a set of annotations (`Snakefile`)
* junnction split-read counts (`Snakefile.study_jxs`)

Gene and exon level counts are produced as full matrices per annotation, per study (if running on an SRA tranche of multiple studies and with multiple annotations, e.g. GencodeV26 & GencodeM23).  These can be used directly with/in recount3.  

Junctions are also produced per-study for recount3.

But junctions are additionally produced for the Snaptron format, i.e. multiple studies are kept in the same Snaptron compilation if doing SRA, also the annotations used to mark junctions as annotated (or novel) is a separate, and typically larger set, than the annotations used in the gene/exon counts.

## Dependencies

The recount-unify Docker image:

https://quay.io/repository/broadsword/recount-unify?tab=tags

Alternatively, you can install/build the following dependencies:

* `Snakemake` is needed for running unifying workflow
* `python3` is used to run the main Snakemake workflow (though `python2` will probably work)
* `pypy` (for python 2.7) is needed for certain steps within the workflow (jx merging and annotating)
* `zstd` is needed for decompressing `recount-pump` outputs
* `pigz` is needed for compressing final outputs
* a wrapper shell script around `pigz` emulating `zcat` called `pcat` in PATH which maps to `pigz --stdout -p2 -d`
* `bgzip` is needed for compressing final, Snaptron-ready outputs (htslib >=1.9)
* `tabix` same as `bgzip` but can be ~any version of htslib
* `sqlite3` same reason as for `bgzip/tabix`
* Python 2.7 with PyLucene 6.5.0 installed (optional, only needed if building Lucene indexes for Snaptron)
* Custom supporting utility binaries need to be built ahead of time in `merge/`, `scripts/`, and `rejoin/`, cd into each and run `make all`
 
Also, the following files/information are needed to run the unifier, but are specific to your project:

* a tab-delimited file of the project study(s) and sample/run ids used in the `recount-pump` run (e.g. `project1.tsv` below)
* compilation ID for your project that doesn't collide with existing recount3/Snaptron2 compilations IDs (`compilation_id`)
* number of samples/runs in the project (`#_samples` below, should be exactly the set of runs/samples which successfully ran through `recount-pump`)

The following are more generic files used across projects sharing the same reference genome (e.g. `hg38`, or `grcm38`) and set of annotations, e.g. `G026,G029,R109,ERCC,SIRV,F006` needed for the gene/exon sums part:

* set of annotation short names used (e.g. `G026,G029,R109,ERCC,SIRV,F006`)
* disjoint exon-to-gene mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed`)
* disjoint exon-to-annotated exon mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed`)
* final gene disjoint-to-annotation mapping file (e.g. `G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed`)
* exon annotation split mapping file `exon_bitmasks.tsv`
* coordinates for exon annotation split mapping file `exon_bitmask_coords.tsv`
* the list of annotated exon intervals used in the `recount-pump` run (e.g. `exons.bed.w_header.gz`)
  this file has to be be exactly the same order as the exon sums that's produced by `bamcount` as part of `recount-pump`
* list of 0's for efficient placeholder for samples w/o any exon sums (`blank_exon_sums`), this should be the same length as the list of annotated exon intervals file (sans header line)
  
The following files are specifically needed only for the junction aggregation part:

* a list of annotated junctions (e.g. `annotated_junctions.tsv.gz`)
* genome reference chromosome sizes file (e.g. `hg38.recount_pump.fa.sizes.tsv`)
* genome reference chromosome FASTA file (e.g. `hg38.recount_pump.fa`), which needs to be exactly the same as what's used in `recount-pump`, available from https://recount-ref.s3.amazonaws.com/hg38/fasta.tar.gz
The chromosome names need to exactly match what's in genome reference chromosome sizes file above, or this file should be regenerated based on the FASTA file.
The genome FASTA and sizes are only needed for generating the dinucleotide splice site motifs for the junctions since STAR does not report the non-canonical splice motifs.

With the exception of the genome FASTA file, the rest of these files can all be downloaded from https://github.com/langmead-lab/monorail-run-data
for recount3 human & mouse projects.  They can also be obtained from https://recount-ref.s3.amazonaws.com/hg38

## Prep
While the Snakemake files do most of the heavy lifting, there are a number of separate steps which are still outside the Snakemake workflows.

This step assumes that all checking/filtering of `recount-pump` output has already been done previously.

First, you need to setup the unifier working subdirectory (of the unifier github checkout):

```mkdir project1
cd project1
python ../sample_ids/assign_compilation_ids.py --accessions-file project1.tsv --compilation-code <compilation_code> > project1.ids.tsv
```

Where `<compilation_code>` (a.k.a. `compilation_id`) needs to checked against the list of existing recount3/Snaptron2 compilations so it doesn't collide, https://github.com/langmead-lab/monorail-run-data/blob/master/recount-unify/sample_ids/all_compilation_codes.tsv

Then, you need to create symlinks from the output of `recount-pump` (assumes the same filesystem):

```/bin/bash -x ../scripts/find_done.sh /path/to/project_recount-pump_output links project_recount-pump_file_prefix```

An example for an human SRA tranche (5) on the IDIES systems:

```/bin/bash -x ../scripts/find_done.sh /full/path/to/langmead/sra_human_v3_5 links sra_human_v3```

Finally, a reference-version specific file named `blank_exon_sums` needs to be in the working directory of the unifier.

This can be symlinked/copied from the directory where the refs/backing files were downloaded to (see Dependencies section).

## Layout of links to recount-pump output

Due to the critical nature of this part for the proper running of the unifier, this get its own section.

The `scripts/find_done.sh` script mentioned in the previous section *should* organize the symlinks to the original recount-pump output directories correctly, however, it's worth checking given that the rest of the Unifier is critically sensitive to how the links are organized.

For example, if you find that you're getting blanks instead of actual integers in the `all.exon_bw_count.pasted.gz` file, it's likely a sign that the input directory hierarchy was not laid out correctly.

Assuming your top level directory for input is called `links`, the expected directory hierarchy for each sequencing run/sample is:

`links/study_loworder/study/run_loworder/run/symlink_to_recount-pump_attempt_directory_for_this_run`

e.g.:

`links/94/SRP019994/83/SRR797083/sra_human_v3_41_in26354_att2`

where `sra_human_v3_41_in26354_att2` is the symlink to the actual recount-pump generated attempt director for run `SRR797083` in study `SRP019994`.

`study_loworder` and `run_loworder` are *always* the last 2 characters of the study and run accessions/IDs respectively.

## Running the Unifier Workflows

### Gene and exon sums
An semi-generic example of the unifier snakemake command for aggregatig gene and exon sums for a human SRA tranche:

```snakemake -j <#_threads> --stats ./stats.json --snakefile ../Snakefile -p --config input=links staging=unified sample_ids_file=project1.ids.tsv annotated_sjs=/path/to/annotated_junctions.tsv.gz existing_sums=/path/to/exons.bed.w_header.gz compilation=sra gene_rejoin_mapping=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.bed exon_rejoin_mapping=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons.bed recount_pump_output=NOT_USED gene_mapping_final=/path/to/G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed exon_bitmasks=/path/to/srav3h.exon_counts.bitmasks.tsv exon_bitmasks_coords=/path/to/srav3h.exon_counts.bitmasks.coords annotation_list=G026,G029,R109,F006,ERCC,SIRV  num_samples=<#samples> num_exons=<#exons>```

`links` is the directory created in the prep step above which contains the symlinks to to all of recount-pump's ourput files.

`unified` is a directory which will be created by the Snakemake workflow which contains intermediate files, once a run's final files have been produced and checked this directory can be removed.

It is imperaptive that the `num_samples=<#samples>` be set to *exactly* the number of samples (sequencing runs) in the tranche/compilation/study being aggregated.
Similarly for `num_exons=<#exons>`.

`num_exons` was 1709834 for human and 710078 for mouse

### Junction counts
An semi-generic example of the unifier snakemake command for aggregating junction counts for a human SRA tranche:

```snakemake -j <#_threads> --stats ./perstudy.jxs.stats.json --snakefile ../Snakefile.study_jxs -p --config input=links staging=unified_jxs annotated_sjs=/path/to/annotated_junctions.hg38.tsv.gz ref_sizes=/path/to/hg38.recount_pump.fa.new_sizes ref_fasta=/path/to/hg38.recount_pump.fa sample_ids_file=project1.ids.tsv  compilation_id=<compilation_id> study_dir=junction_counts_per_study build_sqlitedb=1```

where you need to provide the `compilation_id` used previously to generate the `rail_id`s for this compilation (e.g. `project1.ids.tsv`).

This will generate both the recount3 and Snaptron-ready junction matrices & indices.

For Snaptron you will still need to build the Lucene indices for the sample metadata, please see the last section for details on that which can be applied here.

## QC summary

This is typically run after the unifier has been run on the project/tranche (if you want it to include the intron sums) in the unifier working subdirectory of the project/tranche:

```python3 ../log_qc/parse_logs_for_qc.py --incoming-dir links --sample-mapping ids.tsv --intron-sums intron_counts_summed.tsv > qc.tsv 2> qc.err```

## Super merge of junctions for Snaptron

The super merge of junctions is primarily useful as an input to the srav3h and srav1m Snaptron compilations.
But it can also be used to do other types of cross tranche/compilation merging of junctions (e.g. all human: sra + gtex + tcga).

```merge/super_jx_merge_for_snaptron.sh <prefix> <set_of_annotated_jxs_file.gz> <compilation_id> <num_cpus_for_bgzip>```

example:
```merge/super_jx_merge_for_snaptron.sh srav3_human annotated_junctions.hg38.tsv.gz 11 8```

Will produce 3 final files for Snaptron:

* `junctions.bgz`
* `junctions.bgz.tbi`
* `junctions.sqlite`

`prefix` here is `srav3_human` which acts as the directory prefix for all the tranches, e.g. `srav3_human5` is the 6th human tranche (starts from 0).
The script assumes it will be run in the parent directory where of all tranches' output from recount-unify in appropriately named subdirectories using that prefix (e.g. `srav3_human1`).

You will still need a set of Lucene indexes of your metadata before you have a minimally viable Snaptron compilation.  See the last section in this document for details.

## Generating Lucene indices for sample metadata for Snaptron
You will start with a tab-delimited file of sample metadata.  This file needs to have as the first column the `<rail_id>` generated as part of the unifier run above (typically in `ids.tsv` or `project1.ids.tsv` per tranche, which can be concatented together to form the full samples list for a super merge compilation).

A useable `samples.tsv` file can be as simple as a list of:

`<rail_id>TAB<sample_name>`

or can have hundreds of columns (e.g. TCGA).

Once you have a `samples.tsv` file you're happy with:

`snaptron/deploy/build_lucene_indexes.sh samples.tsv all`

in the same directory where the relevant `junctions.bgz` is located.

NOTE: the Lucene builder script requires that PyLucene 6.5.0 be installed in the python2.7 instance in PATH.

Gene, exon, and base-level Snaptron databases & indexes are not required for a minimally viable Snaptron server and are not covered here.
