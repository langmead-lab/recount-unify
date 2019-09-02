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
#include "ksort.h"
#include "kvec.h"

#include <unordered_map>
template<typename K, typename V>
using hash_map = std::unordered_map<K,V>;

typedef std::vector<std::string> strlist;
typedef std::vector<char*> charlist;
typedef hash_map<std::string, int> strmap;
typedef hash_map<std::string, uint64_t*> countmap;
typedef hash_map<std::string, strmap*> strmap2;

std::string OUTPUT_PREFIX = "exon";

typedef struct {
    char* chrm;
    long start;
    long end;
    char* strand;
} annotation_t;

typedef struct {
    char* chrm;
    long start;
    long end;
    char* strand;
    //raw read count per sample
    uint64_t* counts;
    //original annotation entry
    annotation_t* a;
    const std::string* key;
} annot_heap_t;

//this is the min version
#define heap_lt(a, b) (strcmp((a).chrm, (b).chrm) > 0 || ( strcmp((a).chrm, (b).chrm) == 0 && (a).end > (b).end))
KSORT_INIT(aheap, annot_heap_t, heap_lt)

typedef struct {
    kvec_t(annot_heap_t) heap;
    countmap* seen;
} hs;

//used for annotation-mapping file where we read a 3+ column BED file
static const int CHRM_COL=6;
static const int START_COL=7;
static const int END_COL=8;

static const int GENE_COL=9;
static const int STRAND_COL=11;
static const int COUNTS_START_COL=6;
//determine what type of exon (or if it's an intron) in the sums BED
static const int TYPE_COL=4;
static const int INTRON_TYPE=1;

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

bool DUP_CHECK=true;
bool SKIP_MISSING_EXONS=false;
bool SKIP_PRINTING_GENE_NAMES = false;

//track disjoint exon "keys" (chrm,start,end,name,score,strand) to array indexes
//for actual genes/exons from annotation
//list of actual annotations (annotation_t structs) for gene/exons
typedef hash_map<std::string, annotation_t*> annotation_t_map_t;
typedef hash_map<std::string, strmap*> annotation_map_t;
static const int process_region_line(char* line, const char* delim, annotation_map_t* amap, char*** key_fields, annotation_t_map_t* alist, uint64_t** counts, int last_col, std::string key_column_type, strmap2* glist) {
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
    sprintf(key2_,"%s\t%s\t%s\t%s",(*key_fields)[0],(*key_fields)[1],(*key_fields)[2],(*key_fields)[5]);
    std::string key2(key2_);
   
    //establish map between disjoint exon key and
    //original annotated entity key 
    bool free_key = true;
    if(amap->find(key) == amap->end()) {
        (*amap)[key] = new strmap[1];
        free_key = false;
    }

    //now store the original annotated entity
    bool free_key2 = true;
    auto it = alist->find(key2);
    if(!SKIP_PRINTING_GENE_NAMES && glist->find(key2) == glist->end())
        (*glist)[key2] = new strmap[1];
    if(it == alist->end()) {
        annotation_t* coords = new annotation_t[1];
        coords->chrm = strdup((*key_fields)[0]);
        coords->start = atol((*key_fields)[1]);
        coords->end = atol((*key_fields)[2]);
        coords->strand = strdup((*key_fields)[5]);
        (*alist)[key2] = coords;
        if(!SKIP_PRINTING_GENE_NAMES)
            (*(*glist)[key2])[std::string(strdup(gname))] = 1;
        free_key2 = false;
    }
    else {
        delete key2_;
        key2 = it->first;
    }
    (*(*amap)[key])[key2]=1;
    if(free_key)
        delete key;

    if(line_copy)
		free(line_copy);
	if(line)
		free(line);
    return ret;
}

typedef std::vector<uint64_t*> intlist;
static const int process_counts_line(char* line, const char* delim, annotation_map_t* amap, char*** key_fields, annotation_t_map_t* alist, intlist* counts_list, hs* h, FILE* fout, FILE* ifout, bool dec_start_coord, bool skip_intron_check, strmap* disjoint_exons_seen, strmap2* glist) {
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
    char* key = new char[1024];
	int i = 0;
    int ret = 0;
    int k = 0;
    char* chrm = nullptr;
    long start = -1;
    strmap* annotations = nullptr;
	while(tok != nullptr) {
        if(i <= KEY_FIELD_COL_END) {
            //if intron, print out line, cleanup, and return
            if(!skip_intron_check && i == TYPE_COL && atoi(tok) == INTRON_TYPE) {
                fprintf(ifout, "%s", line);
                counts_list->clear(); 
                delete key;
                if(line_copy)
                    free(line_copy);
                if(line)
                    free(line);
                return 0;
            }
            memcpy((*key_fields)[i],tok,strlen(tok)+1);
            i++;
    		tok = strtok(nullptr, delim);
            continue;
        }
        //make the first (disjoint exon) key
        else if(i == KEY_FIELD_COL_END+1) {
            chrm = (*key_fields)[0];
            start = atol((*key_fields)[1]);
            if(dec_start_coord)
                start--;
            sprintf(key,"%s\t%lu\t%s\t.\t0\t%s",(*key_fields)[0],start,(*key_fields)[2],(*key_fields)[5]);
            //check for duplicates
            if(DUP_CHECK) {
                if(disjoint_exons_seen->find(key) != disjoint_exons_seen->end()) {
                    //bail on this counts line altogether, don't want to overcount 
                    counts_list->clear(); 
                    delete key;
                    fprintf(stderr,"DUP\t%s\n",line);
                    if(line_copy)
                        free(line_copy);
                    if(line)
                        free(line);
                    return -1;
                }
                (*disjoint_exons_seen)[key] = 1;
            }
            if(SKIP_MISSING_EXONS) {
                if(amap->find(key) == amap->end()) {
                    //skip this disjoint exon counts line altogether, rather than error out
                    counts_list->clear(); 
                    delete key;
                    fprintf(stderr,"MISSING_EXON\t%s\n",line);
                    if(line_copy)
                        free(line_copy);
                    if(line)
                        free(line);
                    return -1;
                }
            }
            //get list of annotations (e.g. exons) which overlap this disjoin exon
            annotations = (*amap)[key];
            //fprintf(stderr,"key1:%s\n",key);
            //for each overlapping annotation, we've either added it already
            //or it needs to be added to the heap
            for(auto const& annot_key : *annotations) {
                auto it = h->seen->find(annot_key.first);
                //check if we've already added this annotation to the heap
                if(it == h->seen->end()) {
                    annotation_t* annot_coord = (*alist)[annot_key.first];
                    annot_heap_t aht;
                    aht.chrm = annot_coord->chrm;
                    aht.start = annot_coord->start;
                    aht.end = annot_coord->end;
                    aht.strand = annot_coord->strand;
                    aht.key = &annot_key.first;
                    //TODO replace this with a static buffer of counts arrays, pre-allocated
                    //and a free list
                    aht.counts = new uint64_t[NUM_SAMPLES]();
                    kv_push(annot_heap_t, h->heap, aht);
                    //re-establish heapmin property
                    ks_heapup_aheap(h->heap.n, h->heap.a);    
                    //mark as added to heap
                    (*(h->seen))[annot_key.first] = aht.counts;
                    //need to track counts array for incrementing later in this function
                    counts_list->push_back(aht.counts);
                }
                //added annotation already, but get its counts array
                else
                    counts_list->push_back(it->second);
            }
            //check heap here if the lowest is non-overlapping (by end coordinate) with the current disjoint exon
            int chrm_cmp = strcmp(h->heap.a[0].chrm,chrm);
            while(chrm_cmp != 0 || ( chrm_cmp == 0 && h->heap.a[0].end < start)) {
                //pop top of the heap and print
                annot_heap_t aht = h->heap.a[0];
                //0-based length calc
                long annot_len = aht.end - aht.start;
                if(!SKIP_PRINTING_GENE_NAMES) {
                    bool start = true;
                    for(auto const& kv: *((*glist)[*aht.key])) {
                        if(!start)
                            fprintf(fout,";");
                        start = false;
                        fprintf(fout,"%s",kv.first.c_str());
                    }
                    fprintf(fout,"\t");
                }

                fprintf(fout,"%s\t%lu\t%lu\t%lu\t%s",aht.chrm,aht.start,aht.end,annot_len,aht.strand);
                for(int z = 0; z < NUM_SAMPLES; z++)
                    fprintf(fout,"\t%lu",aht.counts[z]);
                fprintf(fout,"\n");
                h->heap.a[0] = kv_pop(h->heap);
                //maintain heap min property
                ks_heapdown_aheap(0, h->heap.n, h->heap.a);
                //check next chromosome match
                chrm_cmp = strcmp(h->heap.a[0].chrm,chrm);
                //clean up
                h->seen->erase(*aht.key);
                delete aht.counts;
            }
        }
        //reaching here means we're in the counts part of the line
        //k starts at 0 for the counts matrix
        k = i - COUNTS_START_COL;
        //from here for each count, loop through the set of count arrays, update sum at each array position k
        for(auto const& e: *counts_list) {
            //check for sci. notation (typically from R)
            //assumes that the exponent is always two digits, e.g. 1e+05
            int len = strlen(tok);
            if(len > 2 && tok[len-3] == '+')
                e[k] += atof(tok);
            else 
                e[k] += atol(tok);
        }
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
    
static const int read_annotation(FILE* fin, annotation_map_t* amap, annotation_t_map_t* alist, uint64_t*** counts, std::string key_column_type, strmap2* glist) {
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
	    err = process_region_line(strdup(line), "\t", amap, &key_fields, alist, &((*counts)[counts_idx++]), STRAND_COL, key_column_type, glist);
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

void go(std::string annotation_map_file, std::string disjoint_exon_sum_file, std::string key_column_type, bool header, bool dec_start_coord, bool skip_intron_check)
{
    //start with static count matrix of 1024 exons
    uint64_t** counts = nullptr;
    //heap-related based on lh3's BFC's code:
    //https://github.com/lh3/bfc/blob/69ab176e7aac4af482d7d8587e45bfe239d02c96/correct.c
    //kvec_t(annot_heap_t) ah;
    hs h;
    h.heap.n = 0;
    h.heap.m = REGION_BUFFER_SZ;
    h.heap.a = (annot_heap_t*) calloc(sizeof(annot_heap_t),REGION_BUFFER_SZ);
    h.seen = new countmap[1];

    //load mapping from disjoint exons to original annotated exons and genes
    FILE* afp = fopen(annotation_map_file.c_str(), "r");
    annotation_map_t disjoint2annotation; 
    annotation_t_map_t annot2counts;
    strmap2 glist;
    int err = read_annotation(afp, &disjoint2annotation, &annot2counts, &counts, key_column_type, &glist);
    fclose(afp);
   
    //now walk through file of disjoint exon sums 
    FILE* fin = fopen(disjoint_exon_sum_file.c_str(), "r");
    //temporarily holds the distinct fields used for the matching key
    char** key_fields = new char*[KEY_FIELD_COL_END+1];
    for(int i=0;i<=KEY_FIELD_COL_END;i++)
        key_fields[i]=new char[100];
    intlist* counts_list = new intlist[1];

    int last_col = NUM_SAMPLES+GENE_COL-1;
	
    char* line = nullptr;
	size_t length = 0;
	ssize_t bytes_read = getline(&line, &length, fin);
    if(header)
	    bytes_read = getline(&line, &length, fin);
    err = 0;
    //now output annotated sums already calculated
    char* foutname = new char[1024];
    sprintf(foutname,"%s.counts",OUTPUT_PREFIX.c_str());
    FILE* fout = fopen(foutname,"w");
    //save intron sums separately
    char* ifoutname = new char[1024];
    sprintf(ifoutname,"%s.intron_counts",OUTPUT_PREFIX.c_str());
    strmap disjoint_exon_seen;
    FILE* ifout = fopen(ifoutname,"w");
	while(bytes_read != -1) {
        //annotated sums get calculated in this function
	    err = process_counts_line(strdup(line), "\t", &disjoint2annotation, &key_fields, &annot2counts, counts_list, &h, fout, ifout, dec_start_coord, skip_intron_check, &disjoint_exon_seen, &glist);
		bytes_read = getline(&line, &length, fin);
    }
    fclose(ifout);
    for(int i = 0; i < h.heap.n; i++) {
        //heap min property doesn't matter here, just want to get what's left
        annot_heap_t aht = h.heap.a[i];
        //0-based length calc
        long annot_len = aht.end - aht.start;
        if(!SKIP_PRINTING_GENE_NAMES) {
            bool start = true;
            for(auto const& kv: *(glist[*aht.key])) {
                if(!start)
                    fprintf(fout,";");
                start = false;
                fprintf(fout,"%s",kv.first.c_str());
            }
            fprintf(fout,"\t");
        }
        fprintf(fout,"%s\t%lu\t%lu\t%lu\t%s",aht.chrm,aht.start,aht.end,annot_len,aht.strand);
        for(int z = 0; z < NUM_SAMPLES; z++)
            fprintf(fout,"\t%lu",aht.counts[z]);
        fprintf(fout,"\n");
    }
    fclose(fout);
}

int main(int argc, char* argv[]) {
	int o;
	std::string disjoint_exon_sum_file;
	std::string annotation_map_file;
    //this determine which columns in the annotation_map_file
    //are used to form the unique annotation identifier of the annotated entity (exon,gene)
	std::string key_column_type = "exon";
    int num_samples = 0;
    //how large the counts array is, i.e. this is the maximum number of annotated entities
    //expected to be overlapping in the same region at any given time
    int region_buffer_size = 1024;
    bool header = false;
    bool dec_start_coord = false;
    bool skip_intron_check = false;
    bool skip_dup_check = false;
    bool skip_missing_exons = false;
	while((o  = getopt(argc, argv, "a:d:s:k:hcinep:")) != -1) {
		switch(o) 
		{
			case 'a': annotation_map_file = optarg; break;
			case 'd': disjoint_exon_sum_file = optarg; break;
			case 's': num_samples = atoi(optarg); break;
			case 'k': key_column_type = optarg; break;
			case 'b': region_buffer_size = atoi(optarg); break;
            case 'h': header = true; break;
            case 'c': dec_start_coord = true; break;
            case 'i': skip_intron_check = true; break;
            case 'n': skip_dup_check = true; break;
            case 'e': skip_missing_exons = true; break;
			case 'p': OUTPUT_PREFIX = optarg; break;
            case 'g': SKIP_PRINTING_GENE_NAMES = true; break;
		}
	}
	if(num_samples <= 0) {
		std::cerr << "You must give a positive number > 0 for num_samples, currently=" << num_samples << "\n";
		exit(-1);
	}
    NUM_SAMPLES = num_samples;
    REGION_BUFFER_SZ = region_buffer_size;
    DUP_CHECK = !skip_dup_check;
    SKIP_MISSING_EXONS = skip_missing_exons;
	go(annotation_map_file, disjoint_exon_sum_file, key_column_type, header, dec_start_coord, skip_intron_check);
}
