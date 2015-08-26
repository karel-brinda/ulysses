//
//  bloom.c
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//  and Maciek Sykulski
//

#include <string>
#include <boost/regex.hpp>

#include <fstream>
#include <sstream>
//#include <ifstream>    
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "bloom.hpp"
#include "kmer.hpp"

using namespace std;

//Seqid to taxon map
map<string, string> ID_to_taxon_map;

#define _ILLEGAL_NEWICK_CHARS  std::string(":;\\(\\),\\[\\]\t\n\r=")
#define _ILLEGAL_NEWICK_REPLACE std::string("_")
boost::regex re_illegal_newick("["+_ILLEGAL_NEWICK_CHARS+"]");   


KSEQ_INIT(gzFile, gzread)

bloom::bloom():program_version_tag(PROGRAM_VERSION_TAG),nh(0),array(),seed(NULL)
{}

bloom::bloom(coor as_b, unsigned int nh, const char *seedstr):program_version_tag(PROGRAM_VERSION_TAG),nh(0),array(),seed(seedstr) {
    this->init(as_b,nh,seedstr);
}

void bloom::init(coor as_b, unsigned int nh, const char *seedstr) {    
    DF1;
    bloom *bf = this;
    
    seed = seed_t(seedstr);
    
    bf->program_version_tag=PROGRAM_VERSION_TAG;
                
    if (as_b==0){
        bf->nh=0;
        bf->array.clear();
        return ;
    }
    
    assert(nh>0);
    assert(as_b>0);    
    assert(as_b%8==0); // remove this condition?
    
    bf->nh=nh;    
    bf->array.clear();
    bf->array.resize(as_b,false);
    
    DF2;    
}

bloom::~bloom(){
    //DF1;
    //DF2;    
}


int bloom_create(bloom *bf, const char *fn, int both_directions){
    DF1;
    
    gzFile fp;
    kseq_t *seq;
    int l;
           
    int span = bf->seed.span;
    assert(span>0);
    int weight = bf->seed.weight;
    unsigned int nh = bf->nh;
    unsigned int bytes_kmer=compressed_kmer_size(weight);
    
    uchar * compr_kmer = new uchar[bytes_kmer];
    uint64_t * hashes1 = new uint64_t[nh+1];
    
    fp = gzopen(fn, "r");
    assert(fp != Z_NULL);
    seq = kseq_init(fp);
    
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        //fprintf(stderr,"processing: %s\n", seq->name.s);
        
        if(l >= span){
        const uchar *dna=(uchar*)seq->seq.s;        
        //assert(l >= span);
        for (int i=0;i<l-span+1;i++){
            for(int direction=0;direction<=both_directions;direction++){
                
                if(compress_kmer(&dna[i],&bf->seed,bytes_kmer,compr_kmer,direction)!=0){
                    continue;
                }
                
                compute_hashes(compr_kmer, bytes_kmer, hashes1, nh);
                
                for(unsigned int j=0;j<nh;j++){                    
                    bf->array.set(hashes1[j] % bf->array.size(),true);
                }
                
            }
        }
        //printf("seq: %s\n", seq->seq.s);
        }
    }
    //printf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp);
    
    delete [] compr_kmer;
    delete [] hashes1;
    
    DF2;
    return 0;
}

void bloom_print(const bloom *bf, int verbose){
    DF1;

    bitvector::size_type ones = bf->array.count();
    printf("Bloom:\n");
    if(verbose){    
        for(bitvector::size_type i=0;i<bf->array.size();i++){            
            printf("%d",(int)bf->array.test(i));
        }
    }
    printf("\n");
    printf("Ones: %lu\n",ones);
    printf("Density: %f\n",1.0*ones/bf->array.size());

    DF2;
}


int bloom_save(const bloom *bf, const char * filename){
    DF1;
    // make an archive
    std::ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
    oa << (*bf);
    DF2;
    return 0;
}

int bloom_load(bloom *bf, const char * filename){
    DF1;
    // open the archive
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ia(ifs);
    // restore from the archive
    ia >> (*bf);
    DF2;
    return 0;
}


int bloom_shrink(bloom *bf, long factor) {
    DF1;    
    assert(bf->array.size()%(factor*8)==0);  //maybe drop requirement to be exact in bytes?
    bitvector::size_type new_as_b=bf->array.size()/factor;
    bitvector bf_copy;
    bitvector::size_type shift=new_as_b;    
    while(factor>1){
        if (factor%2==1){                        
            factor-=1;
            bf_copy = bf->array;
            bf_copy>>=shift*(factor);
            bf_copy.resize(shift*(factor));
            bf->array.resize(shift*(factor));
            bf->array|=bf_copy;            
        }
        factor/=2;
        bf_copy = bf->array;
        bf_copy>>=shift*(factor);
        bf_copy.resize(shift*(factor));
        bf->array.resize(shift*(factor));
        bf->array|=bf_copy;                
    }    
    DF2;
    return 0;
}


int bloom_or(bloom *bf1, const bloom *bf2){
    DF1;
    
    assert(bf1->array.size()==bf2->array.size());
    
    bf1->array|=bf2->array;      
    
    DF2;
    return 0;
}


int bloom_or(const bloom *bf1, const bloom *bf2, bloom *bf){
    DF1;
    
    assert(bf1->array.size()==bf2->array.size());
    
    bf->nh=bf1->nh;
    bf->seed=bf1->seed;
    
    bf->array=bf1->array;
    bf->array|=bf2->array;    
        
    DF2;
    return 0;
}



int bloom_and(bloom *bf1, const bloom *bf2){
    DF1;
    
    assert(bf1->array.size()==bf2->array.size());
    
    bf1->array&=bf2->array;      
    
    DF2;
    return 0;
}

int bloom_and(const bloom *bf1, const bloom *bf2, bloom *bf){
    DF1;
    
    assert(bf1->array.size()==bf2->array.size());
    
    bf->nh=bf1->nh;
    bf->seed=bf1->seed;
    
    bf->array=bf1->array;
    bf->array&=bf2->array;    
        
    DF2;
    return 0;
}


int bloom_xor(bloom *bf1, const bloom *bf2){
    DF1;
    
    assert(bf1->array.size()==bf2->array.size());
    
    bf1->array^=bf2->array;      
    
    DF2;
    return 0;
}

int bloom_xor(const bloom *bf1, const bloom *bf2, bloom *bf){
    DF1;
    
    assert(bf1->array.size()==bf2->array.size());
    
    bf->nh=bf1->nh;
    bf->seed=bf1->seed;
    
    bf->array=bf1->array;
    bf->array^=bf2->array;      
        
    DF2;
    return 0;
}


bitvector::size_type bloom_ones(const bloom *bf){
    DF1;
    
    bitvector::size_type ones = bf->array.count();
    
    DF2;
    return ones;
}

int bloom_query(const bloom *bf, const uchar *gstr, int direction){
    //DF1;
    
    int span = bf->seed.span;
    assert(span>0);
    int weight = bf->seed.weight;
    coor nh = bf->nh;
    unsigned int bytes_kmer=compressed_kmer_size(weight);
    
    uchar * compr_kmer = new uchar[bytes_kmer];
    uint64_t * hashes1 = new uint64_t[nh+1];
            
    if(compress_kmer(gstr,&bf->seed,bytes_kmer,compr_kmer,direction)!=0){
        return -1;
    }        
    compute_hashes(compr_kmer, bytes_kmer, hashes1, nh);
    
    for(int j=0;j<nh;j++){            
    if (!bf->array.test(hashes1[j] % bf->array.size())){
            return 0;
        }
    }
    
    delete [] compr_kmer;
    delete [] hashes1;
    
    //DF2;        
    return 1;
}



void read_ID_to_taxon_map(const string & ID_to_taxon_map_filename) {
  ifstream map_file(ID_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", ID_to_taxon_map_filename.c_str());
  }
  string line;
  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    string seq_id;
    string taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;    
    ID_to_taxon_map[boost::regex_replace(seq_id,re_illegal_newick,_ILLEGAL_NEWICK_REPLACE)] = taxid;
  }
}



map<string,bloom> * bloom_create_many_blooms(const bloom * initial_bf, const bloom * exclude_bf,
                                             coor as_b, unsigned int nh, 
                                             const char *seedstr, const char *fn, int both_directions){
    DF1;
    
    gzFile fp;
    kseq_t *seq;
    int l;
    
    unsigned int min_nh=nh;
    if (exclude_bf!=nullptr) {
        min_nh=exclude_bf->nh<nh?exclude_bf->nh:nh;
    }
        
    seed_t seed(seedstr);
    int span = seed.span;
    assert(span>0);    
    int weight = seed.weight;
    unsigned int bytes_kmer=compressed_kmer_size(weight);
    
    uchar * compr_kmer = new uchar[bytes_kmer];
    uint64_t * hashes1 = new uint64_t[nh+1];        
    //string taxid;
    
    //Create taxon bloom filter map
    map<string,bloom> * taxon_bloom_map = new map<string,bloom>();
    
    fp = gzopen(fn, "r");
    assert(fp != Z_NULL);
    seq = kseq_init(fp);
    
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        //fprintf(stderr,"processing: %s\n", seq->name.s);
        
        if(l >= span){
        try {
           string & taxid = ID_to_taxon_map.at(boost::regex_replace(std::string(seq->name.s),re_illegal_newick,_ILLEGAL_NEWICK_REPLACE));
           
            //Return already used bloom, or a fresh initialized to 0	            
            bloom & bf = (*taxon_bloom_map)[taxid];	
            
            if (bf.array.size()<=0)  {//This is a new bf, not initialized yet
                if (initial_bf==nullptr)
                    bf.init(as_b, nh, seedstr);
                else
                    bf = *initial_bf;
            }
            
            const uchar *dna=(uchar*)seq->seq.s;
                        
            for (int i=0;i<l-span+1;i++){
                for(int direction=0;direction<=both_directions;direction++){
                    
                    if(compress_kmer(&dna[i],&bf.seed,bytes_kmer,compr_kmer,direction)!=0){
                        continue;
                    }
                    
                    compute_hashes(compr_kmer, bytes_kmer, hashes1, nh);
                    
                    bool do_set = true;
                    if (exclude_bf!=nullptr) {
                       bool all_in = true;                        
                       for(unsigned int j=0;j<min_nh;j++){                    
                         all_in&=exclude_bf->array.test(hashes1[j] % exclude_bf->array.size());
                       }
                       do_set&=!all_in;
                    }
                    if (do_set){
                      for(unsigned int j=0;j<nh;j++){                    
                        bf.array.set(hashes1[j] % bf.array.size(),true);
                      }
                    }
                    
                }
            }
        } 
        catch (const std::out_of_range& oor) {
            fprintf(stderr,"Sequence not mapped to taxid. Omitting.\n");
        }
        }        
        //printf("seq: %s\n", seq->seq.s);
    }
    //printf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp);
    
    delete [] compr_kmer;
    delete [] hashes1;
    
    DF2;
    return taxon_bloom_map;
}
