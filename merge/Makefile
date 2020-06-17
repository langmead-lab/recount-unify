all: sum_counts perbase

sum_counts: sum_counts.c
	gcc -O3 sum_counts.c -o sum_counts

sum_counts-debug: sum_counts.c
	gcc -g sum_counts.c -o sum_counts-debug

perbase: perbase.cpp
	g++ -O3 -std=c++11 perbase.cpp -o perbase

perbase-debug: perbase.cpp
	g++ -g -std=c++11 perbase.cpp -o perbase-debug

clean: 
	rm -f sum_counts sum_counts-debug perbase perbase-debug
