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

#include <stdlib.h>

#include "ulysses/bloom.hpp"
#include "ulysses/kmer.hpp"

using namespace std;

//Seqid to taxon map
map<string, string> ID_to_taxon_map;

#define _ILLEGAL_NEWICK_CHARS  std::string(":;\\(\\),\\[\\]\t\n\r=")
#define _ILLEGAL_NEWICK_REPLACE std::string("_")
boost::regex re_illegal_newick("["+_ILLEGAL_NEWICK_CHARS+"]");   


KSEQ_INIT(gzFile, gzread)

Bloom::Bloom():
    program_version_tag(PROGRAM_VERSION_TAG),nh(0),array(),seed(NULL)
{}

Bloom::Bloom(const Bloom &bf):
    program_version_tag(PROGRAM_VERSION_TAG),nh(bf.nh),array(bf.array),seed(bf.seed)
{
}


Bloom::Bloom(coor as_b, unsigned int nh, const char *seedstr):program_version_tag(PROGRAM_VERSION_TAG),nh(0),array(),seed(seedstr) {
    this->init(as_b,nh,seedstr);
}

void Bloom::init(coor as_b, unsigned int nh, const char *seedstr) {    
    DF1;
    Bloom *bf = this;
    
    seed = Seed(seedstr);
    
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

Bloom::~Bloom(){
    //DF1;
    //DF2;    
}


int Bloom::create(const char *fasta_fn,int both_directions){
    DF1;
    
    gzFile fp;
    kseq_t *seq;
    int l;
           
    int span = this->seed.span;
    assert(span>0);
    int weight = this->seed.weight;
    unsigned int nh = this->nh;
    unsigned int bytes_kmer=compressed_kmer_size(weight);
    
    uchar * compr_kmer = new uchar[bytes_kmer];
    uint64_t * hashes1 = new uint64_t[nh+1];
    
    fp = gzopen(fasta_fn, "r");
    assert(fp != Z_NULL);
    seq = kseq_init(fp);
    
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        //fprintf(stderr,"processing: %s\n", seq->name.s);
        
        if(l >= span){
        const uchar *dna=(uchar*)seq->seq.s;        
        //assert(l >= span);
        for (int i=0;i<l-span+1;i++){
            for(int direction=0;direction<=both_directions;direction++){
                
                if(compress_kmer(&dna[i],&(this->seed),bytes_kmer,compr_kmer,direction)!=0){
                    continue;
                }
                
                compute_hashes(compr_kmer, bytes_kmer, hashes1, nh);
                
                for(unsigned int j=0;j<nh;j++){                    
                    this->array.set(hashes1[j] % this->array.size(),true);
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

void Bloom::print(int verbose) const {
    DF1;

    bitvector::size_type ones = this->array.count();
    printf("Bloom:\n");
    if(verbose){    
        for(bitvector::size_type i=0;i<this->array.size();i++){
            printf("%d",(int)this->array.test(i));
        }
    }
    printf("\n");
    printf("Ones: %lu\n",ones);
    printf("Density: %f\n",1.0*ones/this->array.size());

    DF2;
}


int Bloom::save(const char * filename) const {
    DF1;
    // make an archive
    std::ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
    oa << *this;
    DF2;
    return 0;
}

int Bloom::load(const char * filename){
    DF1;
    // open the archive
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ia(ifs);
    // restore from the archive
    ia >> (*this);
    DF2;
    return 0;
}


int Bloom::shrink(long factor) {
    DF1;    
    assert(this->array.size()%(factor*8)==0);  //maybe drop requirement to be exact in bytes?
    bitvector::size_type new_as_b=this->array.size()/factor;
    bitvector bf_copy;
    bitvector::size_type shift=new_as_b;    
    while(factor>1){
        if (factor%2==1){                        
            factor-=1;
            bf_copy = this->array;
            bf_copy>>=shift*(factor);
            bf_copy.resize(shift*(factor));
            this->array.resize(shift*(factor));
            this->array|=bf_copy;            
        }
        factor/=2;
        bf_copy = this->array;
        bf_copy>>=shift*(factor);
        bf_copy.resize(shift*(factor));
        this->array.resize(shift*(factor));
        this->array|=bf_copy;                
    }    
    DF2;
    return 0;
}

bitvector::size_type Bloom::ones() const {
    DF1;
    
    bitvector::size_type ones = this->array.count();
    
    DF2;
    return ones;
}

int Bloom::query(const uchar *gstr, int direction, unsigned int bytes_compr_kmer) const {
    //DF1;
    
    //uchar * compr_kmer = new uchar[bytes_compr_kmer];
    //uint64_t * hashes1 = new uint64_t[nh+1];
    //Instead of the above: GCC extension: Variable length arrays
    //uchar compr_kmer[bytes_compr_kmer];
    //uint64_t hashes1[nh+1];    
    
    uchar * compr_kmer = (uchar*) alloca(bytes_compr_kmer*sizeof(uchar));
    uint64_t * hashes1 = (uint64_t*) alloca((nh+1)*sizeof(uint64_t));
    
    int span = this->seed.span;
    assert(span>0);    
    coor nh = this->nh;        
            
    if(compress_kmer(gstr,&(this->seed),bytes_compr_kmer,compr_kmer,direction)!=0){
        return -1;
    }        
    compute_hashes(compr_kmer, bytes_compr_kmer, hashes1, nh);
    
    for(int j=0;j<nh;j++){            
        if (!this->array.test(hashes1[j] % this->array.size())){
            return 0;
        }
    }
    
    //DF2;        
    return 1;
}

Bloom& Bloom::operator=(const Bloom &rhs){
    DF1;

    if(rhs.array.size()!=this->array.size()){
        this->array.resize(rhs.array.size());
    }

    this->nh=rhs.nh;
    this->seed=rhs.seed;
    this->array=rhs.array;

    DF2;

    return *this;
}


Bloom& Bloom::operator^=(Bloom const &rhs) {
    DF1;
    
    assert(rhs.array.size()==this->array.size());
    this->array^=rhs.array;      
    
    DF2;

    return *this;
}

Bloom& Bloom::operator|=(Bloom const &rhs) {
    DF1;
    
    assert(rhs.array.size()==this->array.size());
    this->array|=rhs.array;      
    
    DF2;

    return *this;
}

Bloom& Bloom::operator&=(Bloom const &rhs) {
    DF1;
    
    assert(rhs.array.size()==this->array.size());
    this->array&=rhs.array;      
    
    DF2;

    return *this;
}

const Bloom Bloom::operator&(const Bloom &other) const {
    Bloom result = *this;
    result &= other;
    return result;
}

const Bloom Bloom::operator|(const Bloom &other) const {
    Bloom result = *this;
    result |= other;
    return result;
}

const Bloom Bloom::operator^(const Bloom &other) const {
    Bloom result = *this;
    result ^= other;
    return result;
}

bool Bloom::operator==(const Bloom &other) const {
    return this->array == other.array;
}

bool Bloom::operator!=(const Bloom &other) const {
    return !(*this == other);
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



map<string,Bloom> * bloom_create_many_blooms(const Bloom * initial_bf, 
                                             const Bloom * exclude_bf,
                                             const Bloom * include_bf,
                                             const coor as_b, const unsigned int nh, 
                                             const char *seedstr, const char *fn, const int both_directions,
                                             const size_t max_threads,
                                             const size_t work_unit_size){
    DF1;
    
    gzFile fp;
    kseq_t *seq;
    long l;
    
    unsigned int min_nh=nh, min_inc_nh=nh;
    if (exclude_bf!=NULL) {
        min_nh=exclude_bf->nh<nh?exclude_bf->nh:nh;
    }
    if (include_bf!=NULL) {
        min_inc_nh=include_bf->nh<nh?include_bf->nh:nh;
    }
        
    Seed seed(seedstr);
    int span = seed.span;
    assert(span>0);    
    int weight = seed.weight;
    unsigned int bytes_kmer=compressed_kmer_size(weight);
    
    //Create taxon bloom filter map
    map<string,Bloom> * taxon_bloom_map = new map<string,Bloom>();
    
    fp = gzopen(fn, "r");
    assert(fp != Z_NULL);
    seq = kseq_init(fp);
    
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        //fprintf(stderr,"processing: %s\n", seq->name.s);
        
        if(l >= span){
        try {
           string & taxid = ID_to_taxon_map.at(boost::regex_replace(std::string(seq->name.s),re_illegal_newick,_ILLEGAL_NEWICK_REPLACE));
           
            //Return already used bloom, or a fresh initialized to 0	            
            Bloom & bf = (*taxon_bloom_map)[taxid];	
            
            if (bf.array.size()<=0)  {//This is a new bf, not initialized yet
                if (initial_bf==NULL)
                    bf.init(as_b, nh, seedstr);
                else
                    bf = *initial_bf;
            }
            
            const uchar *dna=(uchar*)seq->seq.s;
                        
            long start = 0;
            long end = l-span+1;
            if (!(start<end))
                continue;
            long num_threads = (end-start+work_unit_size-1)/work_unit_size;
            num_threads = num_threads>(long)max_threads?max_threads:num_threads;
            //Enlarge work_unit_size so it's optimal
            long curr_work_unit_size = (end-start+num_threads-1)/num_threads;
            curr_work_unit_size = curr_work_unit_size<(long)work_unit_size ? work_unit_size:curr_work_unit_size;
            #pragma omp parallel num_threads(num_threads)
            {
                uchar * compr_kmer = new uchar[bytes_kmer];
                uint64_t * hashes1 = new uint64_t[nh+1];
                
                long mystart, myend;
                while (true) {
                    #pragma omp critical(get_input_range)
                    {                    
                        mystart = start;
                        myend = end+curr_work_unit_size;
                        myend = myend>end?end:myend;
                        start = myend;                    
                    }// end critical
                    
                    if (!(mystart<myend))
                        break;
                    
                    for (long i=mystart;i<myend;i++){
                        for(int direction=0;direction<=both_directions;direction++){
                            
                            if(compress_kmer(&dna[i],&bf.seed,bytes_kmer,compr_kmer,direction)!=0){
                                continue;
                            }
                            
                            compute_hashes(compr_kmer, bytes_kmer, hashes1, nh);
                            
                            bool do_set = true;
                            if (exclude_bf!=NULL) {
                            bool all_in = true;                        
                            for(unsigned int j=0;j<min_nh;j++){                    
                                all_in&=exclude_bf->array.test(hashes1[j] % exclude_bf->array.size());
                            }
                            do_set&=!all_in;
                            }
                            if (include_bf!=NULL) {
                            bool all_in = true;                        
                            for(unsigned int j=0;j<min_inc_nh;j++){                    
                                all_in&=include_bf->array.test(hashes1[j] % include_bf->array.size());
                            }
                            do_set&=all_in;
                            }
                            if (do_set){
                            for(unsigned int j=0;j<nh;j++){                    
                                bf.array.set(hashes1[j] % bf.array.size(),true);
                            }
                            }
                            
                        }
                    }                
                }  
                delete [] compr_kmer;
                delete [] hashes1;
            } // end parallel
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
    
    DF2;
    return taxon_bloom_map;
}
