//
//  bloom.h
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//  and Maciek Sykulski
//

#ifndef bloom_bloom_h
#define bloom_bloom_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include <time.h>
#include <errno.h>
#include <err.h>
#include <sysexits.h>


#include "boost/dynamic_bitset.hpp"
#include "boost/dynamic_bitset/serialization.hpp"
#include <map>
#include <string>
#include <iostream>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>

#include "utils/MurmurHash3.hpp"
#include "utils/kseq.h"
#include "seed.hpp"
#include "kmer.hpp"
#include "misc.hpp"

#define PROGRAM_VERSION_TAG 2

typedef boost::dynamic_bitset<uint64_t> bitvector;

class bloom {
public:
    bloom();
    bloom(coor as_b, unsigned int nh, const char *seedstr);
    ~bloom();    
    void init(coor as_b, unsigned int nh, const char *seedstr);
    uint32_t program_version_tag;
    uint32_t nh;
    bitvector array;
    seed_t seed;

private:
    friend class boost::serialization::access;
    
    template<class Archive>
    void save(Archive & ar, const unsigned int ) const
    {
        ar & program_version_tag & nh & seed & array;
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int ){
        // invoke serialization of the base class 
        //ar >> boost::serialization::base_object<base_class_of_T>(*this);                
        ar & program_version_tag & nh & seed & array;        
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()    
};

//int bloom_init(bloom *bf, coor as, int nh, const char *seedstr);

int bloom_create(bloom *bf, const char *fasta_fn,int both_directions);

void bloom_print(const bloom *bf, int verbose);

int bloom_load(bloom *bf, const char *fn);

int bloom_save(const bloom *bf, const char *fn);

int bloom_shrink(bloom *bf, long factor);

//int bloom_free(bloom *bf);

int bloom_or(bloom *bf1, const bloom *bf2);
int bloom_or(const bloom *bf1, const bloom *bf2, bloom *bf);

int bloom_xor(bloom *bf1, const bloom *bf2);
int bloom_xor(const bloom *bf1, const bloom *bf2, bloom *bf);

int bloom_and(bloom *bf1, const bloom *bf2);
int bloom_and(const bloom *bf1, const bloom *bf2, bloom *bf);

bitvector::size_type bloom_ones(const bloom *bf);

int bloom_query(const bloom *bf, const uchar *gstr, int dir);

void read_ID_to_taxon_map(const std::string & ID_to_taxon_map_filename);

std::map<std::string,bloom> * bloom_create_many_blooms(const bloom * initial_bf, const bloom * exclude_bf,const bloom * include_bf, coor as_b, unsigned int nh, const char *seedstr, const char *fn, int both_directions);

#endif
