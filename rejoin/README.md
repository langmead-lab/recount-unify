This is for the final step in the recount/snaptron summing pipeline for annotated genes/exons.
`rejoin` is run first twice, to get sums for 1) original annotated exons 2) original annotated genes.
The annotation mapping file that rejoin uses, usually from the monorail-run-data repo but could be custom, must be in coordinate sorted order (typically via `sort -k1,1 -k2,2n -k3,3n`).
Then `rejoin_genes.py` is run to split out genes for the various annotation sources (if there's more than one annotation source present, e.g. Genecode 25 and RefSeq).
