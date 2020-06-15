all: check_for_blanks mmformat

check_for_blanks: check_for_blanks.cpp
	g++ -O3 -std=c++11 check_for_blanks.cpp -o check_for_blanks

check_for_blanks_debug: check_for_blanks.cpp
	g++ -g -std=c++11 check_for_blanks.cpp -o check_for_blanks_debug

mmformat: mmformat.cpp
	g++ -O3 -std=c++11 mmformat.cpp -o mmformat

mmformat-debug: mmformat.cpp
	g++ -g -std=c++11 mmformat.cpp -o mmformat-debug

clean: 
	rm -f check_for_blanks check_for_blanks_debug mmformat mmformat-debug
