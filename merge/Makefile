#all: merge
all: splode

splode: splode.cpp
	#g++ -O3 -std=c++11 splode.cpp -o splode
	g++ -g -std=c++11 splode.cpp -o splode

#merge: merge.cpp
	#g++ -O3 -Wno-write-strings -std=c++11 merge.cpp -o merge
#	g++ -g -Wno-write-strings -std=c++11 merge.cpp -o merge

test: tests/test.bed merge
	cat tests/test.bed | ./merge -c tests/chr1.size -f tests/small_sm.tsv -m -s 5  > test.bed.test
	diff tests/test.bed.expected test.bed.test 

clean: 
	rm merge splode
