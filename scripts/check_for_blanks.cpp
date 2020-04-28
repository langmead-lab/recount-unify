#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <cstdlib>

int NUM_LINES_TO_BUFFER = 10;
int MAX_SIZE_OF_VALUE = 11;

char TAB=9;
char NL=10;

void process_bytes(char* buf, uint32_t bytes_read, bool* prev_newline_tab, uint64_t* line_idx, uint64_t* field_idx)
{
    //TODO: fix col idx reporting, wrong now
    for(uint32_t col_idx = 0; col_idx < bytes_read; col_idx++)
    {
        if(*prev_newline_tab && (buf[col_idx] == NL || buf[col_idx] == TAB))
            fprintf(stdout,"line:%lu,field:%lu\n", *line_idx, *field_idx);
        (*prev_newline_tab) = false;
        if(buf[col_idx] == NL)
        {
            (*line_idx)++;
            (*field_idx) = 0;
            (*prev_newline_tab) = true;
        }
        else if(buf[col_idx] == TAB)
        {
            (*field_idx)++;
            (*prev_newline_tab) = true;
        }
    }
}


int main(int argc, char* argv[]) 
{
    int o;
    uint32_t num_fields = 0;
    uint64_t region_buffer_size = 0;
	while((o  = getopt(argc, argv, "n:b:")) != -1)
    {
		switch(o) 
		{
			case 'n': num_fields = atoi(optarg); break;
			case 'b': region_buffer_size = atol(optarg); break;
		}
    }
    if(num_fields == 0)
    {
        fprintf(stderr,"you must pass the number of fields in each line (assumed to be the same): -n <num_fields>\n");
        exit(-1);
    }
    //setup main buffer ~45 MBs for SRAv1M as default
    if(region_buffer_size == 0)
        region_buffer_size = (NUM_LINES_TO_BUFFER)*num_fields*MAX_SIZE_OF_VALUE;
    //allow some capacity in the buffer for reading to the end of a line
    uint32_t bytes_to_read = region_buffer_size;
    char* buf = new char[region_buffer_size];

    uint64_t line_idx = 0;
    int i = 0;
    uint64_t field_idx = 0;
    bool prev_newline_tab = false;
    uint64_t total_bytes = 0;
    uint32_t bytes_read = fread(buf, 1, bytes_to_read, stdin);
    while(bytes_read == bytes_to_read)
    {
        total_bytes += bytes_read;
        fprintf(stderr,"calling process_bytes running count: %lu\n",total_bytes);
        process_bytes(buf, bytes_read, &prev_newline_tab, &line_idx, &field_idx);
        bytes_read = fread(buf, 1, bytes_to_read, stdin);
    }
    total_bytes += bytes_read;
    fprintf(stderr,"calling process_bytes running count: %lu\n",total_bytes);
    process_bytes(buf, bytes_read, &prev_newline_tab, &line_idx, &field_idx);
    fprintf(stderr,"DONE\n");
}
