## Rejoinning annotation exons/genes

This is for the final step in the recount/snaptron summing pipeline for annotated genes/exons.
`rejoin` is run first twice, to get sums for 1) original annotated exons 2) original annotated genes.
The disjoint exon counts file that rejoin uses must be in coordinate sorted order (typically via `sort -k1,1 -k2,2n -k3,3n`).
Then `rejoin_genes.py` is run to split out genes for the various annotation sources (if there's more than one annotation source present, e.g. Genecode 25 and RefSeq).
`rejoin_genes.sh` calls `rejoin_genes.py` and gzips the output.

## Study/annotation level splitting

* `create_exon_sums_by_study_splits.py` creates the jobs to split the original rejoined exons sums file by `study`

* `split_exons_by_annotation.cpp` (which compiles to `split`) further splits the per-study exons sums into per-stuy, per-annotation exon sums for the final recount3 format

* `split_out_exons.sh` calls `split` for a study and the list of annotations

* `map_exon_sum_rows_to_annotations.py` is used to make the bitmasks files used for in `split`

