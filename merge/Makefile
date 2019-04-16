#all: merge
all: expand

expand: expand.cpp
	g++ -O3 -std=c++11 expand.cpp -o expand
	#g++ -g -std=c++11 expand.cpp -o expand

#merge: merge.cpp
	#g++ -O3 -Wno-write-strings -std=c++11 merge.cpp -o merge
#	g++ -g -Wno-write-strings -std=c++11 merge.cpp -o merge

test: tests/test.bed merge
	cat tests/test.bed | ./merge -c tests/chr1.size -f tests/small_sm.tsv -m -s 5  > test.bed.test
	diff tests/test.bed.expected test.bed.test 

clean: 
	rm merge expand
