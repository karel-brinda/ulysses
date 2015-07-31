//
//  kmer.c
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//
//

#include "kmer.hpp"
#include "seed.hpp"
#include <math.h>
#include "misc.hpp"

unsigned int compressed_kmer_size(unsigned int length){
    return (length +4 -1)/4;
}

int compress_kmer(const uchar *gstr,const seed_t *seed,unsigned int bytes, uchar *compr,int direction) {    
    unsigned int p=0;
    
    for(unsigned int i=0;i<bytes;i++){
        uchar bits=0;
        
        for(unsigned int j=0;j<4;j++){
            bits <<=2;
            
            if(p >= (unsigned int)seed->weight){
                continue;
            }
            
            //printf("%c",gstr[ (seed->shape)[pos] ]);
            switch (gstr[ seed->shape[direction][p++] ]) {
                case 'A':
                case 'a':
                    bits|=0;
                    break;
                case 'C':
                case 'c':
                    bits|=1;
                    break;
                case 'G':
                case 'g':
                    bits|=2;
                    break;
                case 'T':
                case 't':
                    bits|=3;
                    break;
                default:
                    /* not a nucleotide */
                    return 1;
            }
            
            /* complementing (reversing was already done) */
            if (direction==1){
                bits^=3;
            }
        }

        compr[i] = bits;
    }
    
    return 0;
}


int decompress_kmer(const uchar *compr, uchar *gstr, unsigned int bytes) {        
    char de[]={'A','C','G','T'};
    
    unsigned int p=0;
    for (unsigned int i=0;i<bytes;i++){
        for(unsigned int j=0;j<4;j++){
            gstr[p++]=de[3 & (compr[i]>>2*(3-j))];
        }
    }
    
    return 0;
}

int compute_hashes(const uchar *data, const int len, uint32_t *hashes, const int nh) {
    
    for(int j=0;j<nh;j++){
        MurmurHash3_x86_32(data, len, j, &hashes[j]);
    }
    
    return 0;
}

#include <boost/functional/hash.hpp>

int compute_hashes_new(const uchar *data, const int len, uint32_t *hashes, const int nh) {
    
    for(int j=0;j<nh;j++){        
        size_t seed = j;
        for (const uchar *d = data;d<data+len;++d)
            boost::hash_combine(seed,*d);
        hashes[j] = seed^(seed>>32);
    }
    
    return 0;
}
