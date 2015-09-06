//
//  seed.h
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//
//

#ifndef bloom_seed_h
#define bloom_seed_h

#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include <boost/serialization/serialization.hpp>

#define MAX_SEED_SPAN 200

class Seed{
public:
    Seed(const char *seedstr);
    Seed(const Seed &seed);

    char seedstr[MAX_SEED_SPAN];
    int weight;
    int span;
    unsigned char shape[2][MAX_SEED_SPAN];
    
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar,const unsigned int ){
        ar & seedstr & weight & span & shape;
    }
};

#endif
