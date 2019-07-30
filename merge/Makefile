all: rejoin perbase sum_counts

rejoin: rejoin.cpp
	g++ -O3 -std=c++11 rejoin.cpp -o rejoin

rejoin-debug: rejoin.cpp
	g++ -g -std=c++11 rejoin.cpp -o rejoin-debug

perbase: perbase.cpp
	g++ -O3 -std=c++11 perbase.cpp -o perbase

perbase-debug: perbase.cpp
	g++ -g -std=c++11 perbase.cpp -o perbase-debug

sum_counts: sum_counts.c
	gcc -O3 sum_counts.c -o sum_counts

sum_counts-debug: sum_counts.c
	gcc -g sum_counts.c -o sum_counts-debug

test-rejoin: rejoin test/aexon_sums.tsv test/dexon2aexon_map.tsv
	./rejoin -a test/dexon2aexon_map.tsv -d test/aexon_sums.tsv -s 1393 -h -i 2> exon.counts.dups
	sort -k1,4 exon.counts > exon.counts.sorted
	diff test/aexon_sums.rejoined.tsv exon.counts.sorted
	diff test/aexon.counts.dups exon.counts.dups

clean: 
	rm rejoin perbase sum_counts *-debug
