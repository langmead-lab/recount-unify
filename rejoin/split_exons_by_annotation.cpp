#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <string.h>
#include <unistd.h>
#include <cstdlib>

//estimate of the number of columns in the input, used to estimate the line length
//getline will resize as needed
uint32_t EST_NUM_COLS = 50000;

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

int main(int argc, char** argv)
{
	int o;
    //"G026,G029,R109,ERCC,SIRV"
	std::string annotations;
    std::string exon_row_bitmask_file;
    std::string out_prefix;
    bool header = false;
    uint32_t num_rows = 0;
	while((o  = getopt(argc, argv, "a:b:n:p:h")) != -1) {
		switch(o) 
		{
			case 'a': annotations = optarg; break;
			case 'b': exon_row_bitmask_file = optarg; break;
			case 'p': out_prefix = optarg; break;
			case 'n': num_rows = atol(optarg); break;
			case 'h': header = true; break;
        }
    }
	if(annotations.length() == 0 || exon_row_bitmask_file.length() == 0 || num_rows == 0 || out_prefix.length() == 0) {
		std::cerr << "You must pass both: 1) -a \"annotations\" (e.g. G026,G029,R109,ERCC,SIRV) 2) -b exon_sums_row_bitmasks_file and 3) -n <int> the number of rows in the exon sums file 4) -p <prefix> string prefix to be used in the output files (e.g. ERP001942)\n";
		exit(-1);
	}
    char delim = ',';
    std::vector<std::string> aprefixes;
    split_string(annotations, delim, &aprefixes);
    int num_annotations = aprefixes.size();
    //create row bitmask array
    //std::unique_ptr<char[]> bitmasks(new char[num_rows]);
    char* bitmasks = new char[num_rows*num_annotations];
    FILE* fin = fopen(exon_row_bitmask_file.c_str(),"r");
     
	size_t length = -1;
    char* buf = bitmasks;
    uint32_t i = 0;
    ssize_t bytes_read = getline(&buf, &length, fin);
	while(bytes_read != -1)
    {
        i += num_annotations;
        buf = &(bitmasks[i]);
        bytes_read = getline(&buf, &length, fin);
    }
    fprintf(stdout,"num annotations %d\n",num_annotations);
    int j = 0;
    FILE** fps = new FILE*[num_annotations];
    char* fname = new char[100];
    for(j = 0; j < num_annotations; j++)
    {
        sprintf(fname, "%s.%s.tsv", aprefixes[j].c_str(), out_prefix.c_str());
        fps[j] = fopen(fname,"w");
    }
    length = -1;
    i = 0;
    //estimate the line length needed for the buffer
    char* line = new char[EST_NUM_COLS*sizeof(uint32_t)];
    length = EST_NUM_COLS*sizeof(uint32_t);
    bytes_read = getline(&line, &length, stdin);
    if(header)
        bytes_read = getline(&line, &length, stdin);
	while(bytes_read != -1)
    {
        for(j=0; j < num_annotations; j++)
        {
            //check for ASCII "0"
            if(bitmasks[i+j] != 48)
               fprintf(fps[j],"%s",line);
        } 
        i += num_annotations;
        bytes_read = getline(&line, &length, stdin);
    }
    for(j = 0; j < num_annotations; j++)
        fclose(fps[j]);

    //go(num_annotations, aprefixes, num_rows, bitmasks);
    /*char* buf2 = new char[num_annotations+1];
    buf2[num_annotations] = '\0';
    for(i = 0; i < num_rows*num_annotations; i += num_annotations)
    {
        memcpy(buf2, &(bitmasks[i]), num_annotations);
        fprintf(stdout, "%s\n", buf2);
    }
    /*for(int i = 0; i < num_annotations; i++)
        fprintf(stdout, "%s\n", aprefixes[i].c_str());*/

    return 0;
}
