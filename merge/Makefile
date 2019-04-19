all: splode

splode: splode.cpp
	g++ -O3 -std=c++11 splode.cpp -o splode
	#g++ -g -std=c++11 splode.cpp -o splode

#merge: merge.cpp
#	g++ -O3 -Wno-write-strings -std=c++11 merge.cpp -o merge
#	g++ -g -Wno-write-strings -std=c++11 merge.cpp -o merge

#test_merge: tests/test.bed merge
#	cat tests/test.bed | ./merge -c tests/chr1.size -f tests/small_sm.tsv -m -s 5  > test.bed.test
#	diff tests/test.bed.expected test.bed.test 

test: splode
	./splode -a test/dexon2aexon_map.tsv -d test/aexon_sums.tsv -s 1393
	sort -k1,4 exon.counts > exon.counts.sorted
	diff test/aexon_sums.rejoined.tsv exon.counts.sorted

clean: 
	rm merge splode
