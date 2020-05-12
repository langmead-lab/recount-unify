all: rejoin split

rejoin: rejoin.cpp
	g++ -O3 -std=c++11 rejoin.cpp -o rejoin

debug: rejoin.cpp
	g++ -g -std=c++11 rejoin.cpp -o rejoin-debug

split: split_exons_by_annotation.cpp
	g++ -O3 -std=c++11 split_exons_by_annotation.cpp -o split

test: rejoin test/aexon_sums.tsv test/dexon2aexon_map.tsv
	./rejoin -a test/dexon2aexon_map.tsv -d test/aexon_sums.tsv -s 1393 -h -i 2> exon.counts.dups
	sort -k1,4 exon.counts > exon.counts.sorted
	diff test/aexon_sums.rejoined.tsv exon.counts.sorted
	diff test/aexon.counts.dups exon.counts.dups

clean: 
	rm rejoin rejoin-debug split
