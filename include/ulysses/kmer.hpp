//
//  kmer.h
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//
//

#ifndef bloom_kmer_h
#define bloom_kmer_h

#include <stdlib.h>
#include "ulysses/bloom.hpp"
#include "ulysses/seed.hpp"
#include "MurmurHash3.hpp"
#include "ulysses/misc.hpp"

unsigned int compressed_kmer_size(unsigned int length);

int compress_kmer(const uchar *gstr,const Seed *seed, unsigned int bytes, uchar *compr,int direction);

int decompress_kmer(const uchar *compr, uchar *gstr, unsigned int bytes);

int compute_hashes(const uchar *data, const int len, uint64_t *hashes, const int nh);

#endif
