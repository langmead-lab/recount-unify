#include <cstdlib>
#include <fstream>
#include <string>
#include <string.h>
#include <unordered_map>

template<typename K, typename V>
using hash_map = std::unordered_map<K,V>;
typedef hash_map<uint64_t, uint64_t> countmap;

int main(int argc, char** argv)
{
    countmap counts;
    uint64_t num_rows = 0;
    uint64_t sum = 0;
    uint64_t value = 0;
    uint64_t max = 0;
    uint64_t max_count = 0;

    int line_length = 20;
    char* line = new char[line_length];
    size_t length = line_length;
    ssize_t bytes_read = 0;

    while(bytes_read=getline(&line, &length, stdin))
    {
        num_rows++;
        //overwrite newline
        line[bytes_read-1]='\0';
        uint64_t val = strtoull(line, NULL, 0);
        sum += val;
        uint64_t count = ++(counts[val]);
        if(count > max_count)
        {
            max_count = count;
            max = val;
        }
    }
    double avg = sum / num_rows;
    fprintf(stdout,"total_count\t%lu\tsum\t%lu\tavg\t%.3f\tmode\t%lu\tmode_count\t%lu\n", num_rows, sum, avg, max, max_count);
    return 0;
}
