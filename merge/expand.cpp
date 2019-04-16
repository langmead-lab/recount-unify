#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <cstdlib>
#include <assert.h>
#include <vector>
#include <cstdint>
#include <fcntl.h>

#if USE_SKA_FLAT_HASH
    #  include "flat_hash_map.hpp"
    template<typename K, typename V>
    using hash_map = ska::flat_hash_map<K,V>;
#else
#  include <unordered_map>
    template<typename K, typename V>
    using hash_map = std::unordered_map<K,V>;
#endif

//used for annotation-mapping file where we read a 3+ column BED file
static const int CHRM_COL=6;
static const int START_COL=7;
static const int END_COL=8;

static const int GENE_COL=9;

//first 6 columns are the key field for matching between
//the annotation map and the disjoint exon sums
static const int KEY_FIELD_COL_END=5;

//1MB per line should be more than enough for CIO
static const int LINE_BUFFER_LENGTH=1048576;
//this sets the maximum number of overlapping original exons we expect
//in any given region (how many we're prepared to keep in memory until they fall off the heap)
//TODO: replace both of these with user passed in values
static const int MAX_REGION_BUFFER_SZ=500;
//TODO2: figure this out automatically
static const int NUM_SAMPLES=1393;

typedef struct {
    char* name;
    long start;
    long end;
    void* counts;
} annotation_t;

//track disjoint exon "keys" (chrm,start,end,name,score,strand) to array indexes
//for actual genes/exons from annotation
//list of actual annotations (annotation_t structs) for gene/exons
typedef std::vector<char*> charlist;
typedef hash_map<std::string, annotation_t*> annotation_t_map_t;
typedef hash_map<std::string, charlist*> annotation_map_t;
static const int process_region_line(char* line, const char* delim, annotation_map_t* amap, char*** key_fields, annotation_t_map_t* alist) {
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
	int i = 0;

	char* chrm = nullptr;
    long start = -1;
    long end = -1;
	
    char* gname = nullptr;
    
    char* key = new char[1024];

    int ret = 0;
	int last_col = GENE_COL;
	while(tok != nullptr) {
		if(i > last_col)
			break;
        if(i <= KEY_FIELD_COL_END) {
            memcpy((*key_fields)[i],tok,strlen(tok)+1);
            i++;
    		tok = strtok(nullptr, delim);
            continue;
        }
        //make the first (disjoint exon) key
        else if(i == KEY_FIELD_COL_END+1)
            sprintf(key,"%s\t%s\t%s\t%s\t%s\t%s",(*key_fields)[0],(*key_fields)[1],(*key_fields)[2],(*key_fields)[3],(*key_fields)[4],(*key_fields)[5]);
        if(i == CHRM_COL) {
			chrm = strdup(tok);
            memcpy((*key_fields)[i-CHRM_COL],tok,strlen(tok)+1);
        }
        else if(i == START_COL) {
			start = atol(tok);
            memcpy((*key_fields)[i-CHRM_COL],tok,strlen(tok)+1);
        }
        else if(i == END_COL) {
			end = atol(tok);
            memcpy((*key_fields)[i-CHRM_COL],tok,strlen(tok)+1);
        }
        else if(i == GENE_COL) {
			gname = strdup(tok);
            memcpy((*key_fields)[i-CHRM_COL],tok,strlen(tok)+1);
        }
		i++;
		tok = strtok(nullptr, delim);
	}
    //make 2nd key (original annotated exon)
    char* key2 = new char[1024];
    sprintf(key2,"%s\t%s\t%s\t%s",(*key_fields)[0],(*key_fields)[1],(*key_fields)[2],(*key_fields)[3]);
   
    //establish map between disjoint exon key and
    //original annotated exon key 
    bool free_key = true;
    if(amap->find(key) == amap->end()) {
        (*amap)[key] = new charlist[1];
        free_key = false;
    }
    (*amap)[key]->push_back(key2);
    if(free_key)
        delete key;

    //now store the original annotated exon
    free_key = true;
    if(alist->find(key2) == alist->end()) {
        annotation_t* coords = new annotation_t[1];
        coords->start = start;
        coords->end = end;
        coords->name = gname;
        coords->counts = nullptr;
        (*alist)[key2] = coords;
        free_key = false;
    }
    if(free_key)
        delete key2;

    //fprintf(stderr,"key1 %s key2 %s\n",key,key2);

    if(line_copy)
		free(line_copy);
	if(line)
		free(line);
    return ret;
}
    
static const int read_annotation(FILE* fin, annotation_map_t* amap) {
    //temporarily holds the distinct fields used for the matching key
    char** key_fields = new char*[KEY_FIELD_COL_END+1];
    for(int i=0;i<=KEY_FIELD_COL_END;i++)
        key_fields[i]=new char[100];
	char* line = new char[LINE_BUFFER_LENGTH];
	size_t length = LINE_BUFFER_LENGTH;
	assert(fin);
	ssize_t bytes_read = getline(&line, &length, fin);
	int err = 0;
    annotation_t_map_t alist;
	while(bytes_read != -1) {
	    err = process_region_line(strdup(line), "\t", amap, &key_fields, &alist);
        assert(err==0);
		bytes_read = getline(&line, &length, fin);
    }
    delete[] key_fields;
	fprintf(stderr,"building annotation set done, disjoint2annot map size: %lu, original annotation map size: %lu\n", amap->size(), alist.size());
    /*for(auto const& kv : *amap)
        fprintf(stderr,"%s %lu\n", kv.first.c_str(), (kv.second)->size());*/
    return err;
}

void go(std::string annotation_map_file, std::string disjoint_exon_sum_file)
{
    FILE* afp = fopen(annotation_map_file.c_str(), "r");
    annotation_map_t annotations; 
    int err = read_annotation(afp, &annotations);
    fclose(afp);
}

int main(int argc, char* argv[])
{
	int o;
	std::string disjoint_exon_sum_file;
	std::string annotation_map_file;
	while((o  = getopt(argc, argv, "a:d:")) != -1) 
	{
		switch(o) 
		{
			case 'a': annotation_map_file = optarg; break;
			case 'd': disjoint_exon_sum_file = optarg; break;
		}
	}
	go(annotation_map_file, disjoint_exon_sum_file);
}
