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
static const int STRAND_COL=11;
static const int COUNTS_START_COL=6;

//first 6 columns are the key field for matching between
//the annotation map and the disjoint exon sums
static const int KEY_FIELD_COL_END=5;

//1MB per line should be more than enough for CIO
static const int LINE_BUFFER_LENGTH=1048576;
//this sets the maximum number of overlapping original exons we expect
//in any given region (how many we're prepared to keep in memory until they fall off the heap)
//both of these are overwritten with passed in values
static int REGION_BUFFER_SZ=0;
//TODO2: figure this out automatically
static int NUM_SAMPLES=0;

typedef std::vector<std::string> strlist;
typedef std::vector<char*> charlist;
typedef hash_map<std::string, int> charmap;
typedef hash_map<std::string, uint32_t*> countmap;

typedef struct {
    //gene names that have this exon
    charmap* names;
    //raw read count per sample
    uint32_t* counts;
} annotation_t;

//track disjoint exon "keys" (chrm,start,end,name,score,strand) to array indexes
//for actual genes/exons from annotation
//list of actual annotations (annotation_t structs) for gene/exons
typedef hash_map<std::string, annotation_t*> annotation_t_map_t;
typedef hash_map<std::string, charmap*> annotation_map_t;
static const int process_region_line(char* line, const char* delim, annotation_map_t* amap, char*** key_fields, annotation_t_map_t* alist, uint32_t** counts, int last_col, std::string key_column_type) {
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
	int i = 0;

	char* chrm = nullptr;
    long start = -1;
    long end = -1;
	
    char* gname = nullptr;
    char* strand = nullptr;
    
    char* key = new char[1024];

    int ret = 0;
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
        else if(i == STRAND_COL) {
            strand = strdup(tok);
            memcpy((*key_fields)[i-CHRM_COL],tok,strlen(tok)+1);
        }

		i++;
		tok = strtok(nullptr, delim);
	}
    //make 2nd key (original annotated exon/gene/other)
    char* key2_ = new char[1024];
    if(strcmp(key_column_type.c_str(),"exon") == 0)
        sprintf(key2_,"%s\t%s\t%s\t%s",(*key_fields)[0],(*key_fields)[1],(*key_fields)[2],(*key_fields)[5]);
    //gene
    else
        sprintf(key2_,"%s",(*key_fields)[3]);
    std::string key2(key2_);
   
    //establish map between disjoint exon key and
    //original annotated entity key 
    bool free_key = true;
    if(amap->find(key) == amap->end()) {
        (*amap)[key] = new charmap[1];
        free_key = false;
    }

    //now store the original annotated entity
    bool free_key2 = true;
    auto it = alist->find(key2);
    if(it == alist->end()) {
        annotation_t* coords = new annotation_t[1];
        coords->names = new charmap[1];
        coords->counts = *counts;
        (*alist)[key2] = coords;
        free_key2 = false;
    }
    else {
        delete key2_;
        key2 = it->first;
    }
    (*(*amap)[key])[key2]=1;
    if(free_key)
        delete key;
    //add gene name to exon's gene map
    if(strcmp(key_column_type.c_str(),"exon") == 0)
        (*(*alist)[key2]->names)[gname]=1;

    if(line_copy)
		free(line_copy);
	if(line)
		free(line);
    return ret;
}

typedef std::vector<uint32_t*> intlist;
static const int process_counts_line(char* line, const char* delim, annotation_map_t* amap, char*** key_fields, annotation_t_map_t* alist, intlist* counts_list) {
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
    char* key = new char[1024];
	int i = 0;
    int ret = 0;
    int k = 0;
    charmap* annotations = nullptr;
	while(tok != nullptr) {
        if(i <= KEY_FIELD_COL_END) {
            memcpy((*key_fields)[i],tok,strlen(tok)+1);
            i++;
    		tok = strtok(nullptr, delim);
            continue;
        }
        //make the first (disjoint exon) key
        else if(i == KEY_FIELD_COL_END+1) {
            sprintf(key,"%s\t%s\t%s\t%s\t%s\t%s",(*key_fields)[0],(*key_fields)[1],(*key_fields)[2],(*key_fields)[3],(*key_fields)[4],(*key_fields)[5]);
            annotations = (*amap)[key];
            for(auto const& kv: *annotations) {
                uint32_t* c = (*alist)[kv.first]->counts;
                counts_list->push_back(c);
            }
        }
        //reaching here means we're in the counts part of the line
        //k starts at 0 for the counts matrix
        k = i - COUNTS_START_COL;
        //from here for each count, loop through the set of count arrays, update sum at each array position k
        for(auto const& e: *counts_list)
            e[k]+=atol(tok);
		i++;
		tok = strtok(nullptr, delim);
    }
    //need this for each line independently, so make sure to reset 
    counts_list->clear(); 
    delete key;
    if(line_copy)
		free(line_copy);
	if(line)
		free(line);
    return ret;
}
    
static const int read_annotation(FILE* fin, annotation_map_t* amap, annotation_t_map_t* alist, uint32_t*** counts, std::string key_column_type) {
    //temporarily holds the distinct fields used for the matching key
    char** key_fields = new char*[KEY_FIELD_COL_END+1];
    for(int i=0;i<=KEY_FIELD_COL_END;i++)
        key_fields[i]=new char[100];
	char* line = new char[LINE_BUFFER_LENGTH];
	size_t length = LINE_BUFFER_LENGTH;
	assert(fin);
	ssize_t bytes_read = getline(&line, &length, fin);
	int err = 0;
    int counts_idx = 0;
	while(bytes_read != -1) {
	    err = process_region_line(strdup(line), "\t", amap, &key_fields, alist, &((*counts)[counts_idx++]), STRAND_COL, key_column_type);
        assert(err==0);
		bytes_read = getline(&line, &length, fin);
    }
    delete[] key_fields;
	fprintf(stderr,"building annotation set done, disjoint2annot map size: %lu, original annotation map size: %lu\n", amap->size(), alist->size());
    /*for(auto const& kv : *amap) {
        fprintf(stderr,"%s %lu\n", kv.first.c_str(), (kv.second)->size());
        for(auto const& kv2 : *kv.second)
            fprintf(stderr,"\t%s\n",kv2.first.c_str());
    }*/
    return err;
}

void go(std::string annotation_map_file, std::string disjoint_exon_sum_file, std::string key_column_type)
{
    //start with static count matrix of 1024 exons
    uint32_t** counts = new uint32_t*[REGION_BUFFER_SZ];
    for(int i = 0; i < REGION_BUFFER_SZ; i++)
        counts[i] = new uint32_t[NUM_SAMPLES]();

    //load mapping from disjoint exons to original annotated exons and genes
    FILE* afp = fopen(annotation_map_file.c_str(), "r");
    annotation_map_t disjoint2annotation; 
    annotation_t_map_t annot2counts;
    int err = read_annotation(afp, &disjoint2annotation, &annot2counts, &counts, key_column_type);
    fclose(afp);
   
    //now walk through file of disjoint exon sums 
    FILE* fin = fopen(disjoint_exon_sum_file.c_str(), "r");
    //temporarily holds the distinct fields used for the matching key
    char** key_fields = new char*[KEY_FIELD_COL_END+1];
    for(int i=0;i<=KEY_FIELD_COL_END;i++)
        key_fields[i]=new char[100];

	char* line = nullptr;
	size_t length = 0;
	ssize_t bytes_read = getline(&line, &length, fin);
    err = 0;
    int counts_idx = 0;
    int last_col = NUM_SAMPLES+GENE_COL-1;
    uint32_t* counts_temp = new uint32_t[NUM_SAMPLES]();
    intlist* counts_list = new intlist[1];
	while(bytes_read != -1) {
        //annotated exon sums get calculated in this function
	    err = process_counts_line(strdup(line), "\t", &disjoint2annotation, &key_fields, &annot2counts, counts_list);
        assert(err==0);
		bytes_read = getline(&line, &length, fin);
    }
    //now output annotated sums alread calculated
    char* foutname = new char[1024];
    sprintf(foutname,"%s.counts",key_column_type.c_str());
    FILE* fout = fopen(foutname,"w");
    for(auto const& annot : annot2counts) {
        uint32_t* c = annot2counts[annot.first]->counts;
        fprintf(fout,"%s",annot.first.c_str());
        for(int i=0; i < NUM_SAMPLES; i++)
            fprintf(fout,"\t%u",c[i]);
        fprintf(fout,"\n");
    }
    fclose(fout);
    delete[] key_fields;
}

int main(int argc, char* argv[]) {
	int o;
	std::string disjoint_exon_sum_file;
	std::string annotation_map_file;
    //this determine which columns in the annotation_map_file
    //are used to form the unique annotation identifier of the annotated entity (exon,gene)
	std::string key_column_type = "exon";
    int num_samples = 0;
    //how large the counts array is, i.e. this is the maximum number of annated entities
    //expected to be overlapping in the same region at any given time
    int region_buffer_size = 1024;
	while((o  = getopt(argc, argv, "a:d:s:k:")) != -1) {
		switch(o) 
		{
			case 'a': annotation_map_file = optarg; break;
			case 'd': disjoint_exon_sum_file = optarg; break;
			case 's': num_samples = atoi(optarg); break;
			case 'k': key_column_type = optarg; break;
			case 'b': region_buffer_size = atoi(optarg); break;
		}
	}
	if(num_samples <= 0) {
		std::cerr << "You must give a positive number > 0 for num_samples, currently=" << num_samples << "\n";
		exit(-1);
	}
    NUM_SAMPLES = num_samples;
    REGION_BUFFER_SZ = region_buffer_size;
	go(annotation_map_file, disjoint_exon_sum_file, key_column_type);
}
