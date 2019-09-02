This is for the final step in the recount/snaptron summing pipeline for annotated genes/exons.
`rejoin` is run first twice, to get sums for 1) original annotated exons 2) original annotated genes.
Then `rejoin_genes.py` is run to split out genes for the various annotation sources (if there's more than one annotation source present, e.g. Genecode 25 and RefSeq).
