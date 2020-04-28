all: check_for_blanks

check_for_blanks: check_for_blanks.cpp
	g++ -O3 -std=c++11 check_for_blanks.cpp -o check_for_blanks

debug: check_for_blanks.cpp
	g++ -g -std=c++11 check_for_blanks.cpp -o check_for_blanks.debug

clean: 
	rm -f check_for_blanks check_for_blanks.debug
