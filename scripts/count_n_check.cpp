#include <iostream>
#include <string>
#include <vector>
#include <string.h>
#include <unistd.h>

//check output for mmformat, specifically 1) number of cols (e.g. wc -l) *and* 2) range of col index


int main(int argc, char** argv)
{
    int o;
    uint64_t max_col_id = 0;
    while((o  = getopt(argc, argv, "m:")) != -1) 
    {
        switch(o) 
        {
            case 'm': max_col_id = strtoull(optarg,nullptr,0); break;
        }
    }
    if(max_col_id == 0) 
    {
        std::cerr << "You must pass the maximum column index you expect\n";
        exit(-1);
    }
  
    size_t length = -1;
    uint32_t read_size = 1024;
    char* buf = new char[read_size];
    //expects sample_id format to be: ",sample_id:"
    ssize_t bytes_read = getline(&buf, &length, stdin);
    uint64_t num_samples = 0;
    uint64_t num_max_col_ids = 0;
    uint64_t col_id = 0;
    while(bytes_read != -1)
    {
        num_samples++;
        //overwrite newline
        buf[bytes_read-1] = '\0';
        col_id = strtoull(buf,nullptr,0);
        if(col_id > max_col_id) {
            fprintf(stderr,"ERROR\tcolumn index > max column index: %lu > %lu\n",col_id,max_col_id);
            exit(-1);
        }
        if(col_id == max_col_id)
            num_max_col_ids++;
        bytes_read = getline(&buf, &length, stdin);
    }
    fprintf(stdout,"%lu\n",num_max_col_ids);
    fprintf(stdout,"%lu\n",num_samples);
}
