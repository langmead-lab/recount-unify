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

//max number of bytes for a region spec., e.g. chr1|1|10
int SIZE_OF_COORDINATES = 45;

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
    std::string exon_row_coords_file;
    std::string out_prefix;
    bool header = false;
    uint32_t num_rows = 0;
    int num_digits_per_annotation = 3;
	while((o  = getopt(argc, argv, "a:b:n:p:c:d:h")) != -1) {
		switch(o) 
		{
			case 'a': annotations = optarg; break;
			case 'b': exon_row_bitmask_file = optarg; break;
			case 'c': exon_row_coords_file = optarg; break;
			case 'p': out_prefix = optarg; break;
			case 'n': num_rows = atol(optarg); break;
			case 'd': num_digits_per_annotation = atoi(optarg); break;
			case 'h': header = true; break;
        }
    }
	if(annotations.length() == 0 || exon_row_bitmask_file.length() == 0 || exon_row_coords_file.length() == 0 || num_rows == 0 || out_prefix.length() == 0) {
		std::cerr << "You must pass both: 1) -a \"annotations\" (e.g. G026,G029,R109,F006,ERCC,SIRV) 2) -b exon_sums_row_bitmasks_file 3) -c exon_sums_row_coords_file and 4) -n <int> the number of rows in the exon sums file 4) -p <prefix> string prefix to be used in the output files (e.g. ERP001942)\n";
		exit(-1);
	}
    char delim = ',';
    std::vector<std::string> aprefixes;
    split_string(annotations, delim, &aprefixes);
    int num_annotations = aprefixes.size();
    //create row bitmask array
    //std::unique_ptr<char[]> bitmasks(new char[num_rows]);
    
    //***read in bitmasks
    char* bitmasks = new char[num_rows*num_annotations*num_digits_per_annotation];
    FILE* fin = fopen(exon_row_bitmask_file.c_str(),"r");
     
	size_t length = -1;
    char* buf = bitmasks;
    uint32_t i = 0;
    ssize_t bytes_read = getline(&buf, &length, fin);
    int offset = num_annotations*num_digits_per_annotation;
	while(bytes_read != -1)
    {
        i += offset;
        buf = &(bitmasks[i]);
        bytes_read = getline(&buf, &length, fin);
    }
    fclose(fin);
    fprintf(stdout,"num annotations %d\n",num_annotations);
    
    //***read in exon row coordinates (recount_exon_id)
    //expects a single tab delimited field 
    //where the field may be internally delimited by something else (e.g. "|")
    fin = fopen(exon_row_coords_file.c_str(),"r");
    char** exon_row_coords_ptrs = new char*[num_rows];
    char* exon_row_coords = new char[num_rows*SIZE_OF_COORDINATES];

    length = -1;
    buf = exon_row_coords;
    uint32_t exon_row_coords_ptrs_idx = 0;
    uint64_t exon_row_coords_idx = 0;
    bytes_read = getline(&buf, &length, fin);
	while(bytes_read != -1)
    {
        //set ptr in num_rows sized array
        exon_row_coords_ptrs[exon_row_coords_ptrs_idx++] = &(exon_row_coords[exon_row_coords_idx]);
        exon_row_coords_idx += bytes_read;
        //overwrite newline null
        exon_row_coords[exon_row_coords_idx-1] = '\0';

        //get next row
        buf = &(exon_row_coords[exon_row_coords_idx]);
        bytes_read = getline(&buf, &length, fin);
    }
    fclose(fin); 
    fprintf(stdout,"num coordinate bytes read %lu, num ptrs %d\n",exon_row_coords_idx, exon_row_coords_ptrs_idx);

    //***open output files to write rows into
    int j = 0;
    FILE** fps = new FILE*[num_annotations];
    char* fname = new char[100];
    for(j = 0; j < num_annotations; j++)
    {
        sprintf(fname, "%s.%s.tsv", aprefixes[j].c_str(), out_prefix.c_str());
        fps[j] = fopen(fname,"w");
    }

    //***setup vars for actual splitting
    length = -1;
    i = 0;
    exon_row_coords_ptrs_idx = 0;
    //estimate the line length needed for the buffer
    char* line = new char[EST_NUM_COLS*sizeof(uint32_t)];
    length = EST_NUM_COLS*sizeof(uint32_t);
    bytes_read = getline(&line, &length, stdin);
    if(header)
        bytes_read = getline(&line, &length, stdin);
    
    //***now do actual line splitting
    int fp_idx = 0;
    int val = 0;
    int z = 0;
    char* cval = new char[num_digits_per_annotation+1];
    cval[num_digits_per_annotation] = '\0';
    //temp buf for tokenizing each mask
	while(bytes_read != -1)
    {
        fp_idx = 0;
        for(j=0; j < offset; j+=num_digits_per_annotation)
        {
            memcpy(cval, &(bitmasks[i+j]), num_digits_per_annotation);
            val = atoi(cval);
            if(val != 0)
            {
                //print the same line "val" number of times
                //this is due to duplicate exons in annotations shared by different genes
                for(z = 0; z < val; z++)
                    fprintf(fps[fp_idx], "%s\t%s", exon_row_coords_ptrs[exon_row_coords_ptrs_idx], line);
            }
            fp_idx++;
        } 
        i += offset;
        exon_row_coords_ptrs_idx++;
        bytes_read = getline(&line, &length, stdin);
    }
    for(j = 0; j < num_annotations; j++)
        fclose(fps[j]);

    return 0;
}
