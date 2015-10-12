//
//  main.h
//  bloom
//
//  Created by Karel Brinda on 30.04.15.
//  and Maciek Sykulski
//

#ifndef bloom_main_h
#define bloom_main_h

#define DEFAULT_SEED "###-#-##-#-##--##--###-#-#####"

#include <stdio.h>
#include <math.h>

#include "ulysses/bloom.hpp"
#include "ulysses/seed.hpp"
#include "MurmurHash3.hpp"
#include "kseq.h"

int main_create(int argc,char** argv);

int main_create_many(int argc,char** argv);

int main_stats(int argc,char** argv);

int main_mode(int argc,char** argv);

int main_bitwise(int argc,char** argv);

int main_shrink(int argc,char** argv);

int main_dump(int argc,char** argv);

int main_hamming(int argc,char** argv);

int main_symmdiffmat(int argc,char** argv);

int main_query(int argc,char** argv);

int main_query_and_split(int argc,char** argv);

int usage();

#endif
