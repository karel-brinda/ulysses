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

#include "MurmurHash3.hpp"
#include "kseq.h"

#include "ulysses/seed.hpp"
#include "ulysses/kmer.hpp"
#include "ulysses/misc.hpp"

#define PROGRAM_VERSION_TAG 2

typedef boost::dynamic_bitset<uint64_t> bitvector;

class Bloom {
public:
    Bloom();
    Bloom(const Bloom &bf);
    Bloom(coor as_b, unsigned int nh, const char *seedstr);
    ~Bloom();

    void init(coor as_b, unsigned int nh, const char *seedstr);

    int create(const char *fasta_fn,int both_directions);
    void print(int verbose) const;
    int load(const char *fn);
    int save(const char *fn) const;
    int shrink(long factor);
    int query(const uchar *gstr, int dir, unsigned int bytes_kmer) const;
    bitvector::size_type ones() const;

    Bloom& operator=(const Bloom &rhs);

    Bloom& operator^=(Bloom const &rhs);
    Bloom& operator&=(Bloom const &rhs);
    Bloom& operator|=(Bloom const &rhs);

    const Bloom operator^(const Bloom &other) const;
    const Bloom operator&(const Bloom &other) const;
    const Bloom operator|(const Bloom &other) const;

    bool operator==(const Bloom &other) const;
    bool operator!=(const Bloom &other) const;

    //const Bloom operator!() const;

    uint32_t program_version_tag;
    uint32_t nh;
    bitvector array;
    Seed seed;

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

void read_ID_to_taxon_map(const std::string & ID_to_taxon_map_filename);

std::map<std::string,Bloom>* bloom_create_many_blooms(const Bloom * initial_bf, const Bloom * exclude_bf,const Bloom * include_bf, coor as_b, unsigned int nh, const char *seedstr, const char *fn, int both_directions);

#endif
