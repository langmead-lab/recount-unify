dir=$(dirname $0)
#typically ids.tsv but might be different
ids_mapping=$1

zcat all.exon_bw_count.pasted.gz | head -1 | cut -f 7- > all.exon_bw_count.pasted.gz.samples_header

cat all.gene_counts.rejoined.tsv | tail -n+2 | pypy ${dir}/rejoin_genes.py ../G029.G026.R109.F006.20190220.gtf.disjoint2exons2genes.rejoin_genes.bed gene all.exon_bw_count.pasted.gz.samples_header G026,G029,R109,F006,ERCC,SIRV $ids_mapping G026

for t in G026 G029 R109 ERCC SIRV; do
    gzip ${t}.gene.sums.tsv &
done
gzip F006.gene.sums.tsv 
