all: check_for_blanks mmformat count_n_check

check_for_blanks: check_for_blanks.cpp
	g++ -O3 -std=c++11 check_for_blanks.cpp -o check_for_blanks

check_for_blanks_debug: check_for_blanks.cpp
	g++ -g -std=c++11 check_for_blanks.cpp -o check_for_blanks_debug

mmformat: mmformat.cpp
	g++ -O3 -std=c++11 mmformat.cpp -o mmformat

mmformat-debug: mmformat.cpp
	g++ -g -std=c++11 mmformat.cpp -o mmformat-debug

count_n_check: count_n_check.cpp
	g++ -O3 -std=c++11 count_n_check.cpp -o count_n_check

clean: 
	rm -f check_for_blanks check_for_blanks_debug mmformat mmformat-debug count_n_check
