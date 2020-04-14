#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <unistd.h>
#include <cstdlib>

//multiple sources for this kind of tokenization, one which was useful was:
//https://yunmingzhang.wordpress.com/2015/07/14/how-to-read-file-line-by-lien-and-split-a-string-in-c/
void split_string(std::string line, char delim, std::vector<std::string>* tokens)
{
	tokens->clear();
	std::stringstream ss(line);
	std::string token;
	while(getline(ss, token, delim))
	{
		tokens->push_back(token);
	}
}


//input line from rejoined exon counts file:
//ERCC-00002      ERCC-00002      0       1061    1061    + count1 count2 count3...

int main(char** argv, int argc)
{
	int o;
    //"G026,G029,R109,ERCC,SIRV"
	std::string annotations;
    std::string exon_row_bitmask_file;
	while((o  = getopt(argc, argv, "a:b:")) != -1) {
		switch(o) 
		{
			case 'a': annotations = optarg; break;
			case 'b': exon_row_bitmask_file = optarg; break;
        }
    }
	if(annotations.length() == 0 || exon_row_bitmask_file.length() == 0) {
		std::cerr << "You must pass both: 1) -a \"annotations\" (e.g. G026,G029,R109,ERCC,SIRV) and 2) -b exon_sums_row_bitmasks_file\n";
		exit(-1);
	}
    char delim = ',';
    std::vector<std::string> aprefixes;
    split_string(annotations, delim, &aprefixes);
    int num_annotations = aprefixes.length()
    //go(annotations, exon_row_bitmasks);
    return 0;
}
