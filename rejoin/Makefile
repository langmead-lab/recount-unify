all: rejoin split_exons

rejoin: rejoin.cpp
	g++ -O3 -std=c++11 rejoin.cpp -o rejoin

debug: rejoin.cpp split_exons_by_annotation.cpp
	g++ -g -std=c++11 rejoin.cpp -o rejoin-debug
	g++ -g -std=c++11 split_exons_by_annotation.cpp -o split_exons-debug

split_exons: split_exons_by_annotation.cpp
	g++ -O3 -std=c++11 split_exons_by_annotation.cpp -o split_exons

test: rejoin test/aexon_sums.tsv test/dexon2aexon_map.tsv test/aexon_sums.rejoined.tsv test/aexon.counts.dups
	./rejoin -a test/dexon2aexon_map.tsv -d test/aexon_sums.tsv -s 1393 -h -i 2> exon.counts.dups
	sort -k1,4 exon.counts > exon.counts.sorted
	/bin/bash -c "diff <(cut -f 2- test/aexon_sums.rejoined.tsv) <(cut -f 2- exon.counts.sorted)"
	diff test/aexon.counts.dups exon.counts.dups

clean: 
	rm -f rejoin rejoin-debug split_exons split_exons-debug
