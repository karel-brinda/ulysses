//
//  misc.h
//  ulysses
//
//  Created by Karel Brinda on 01.05.15.
//
//

#ifndef ulysses_misc_h
#define ulysses_misc_h

#define DEBUG 1


#include <stdint.h>

typedef int64_t coor;
typedef unsigned char uchar;

#if DEBUG
#define DF1 \
fprintf(stderr,"\x1B[32m"); \
clock_t t = clock(); \
fprintf(stderr,"%s started\n",__PRETTY_FUNCTION__); \
fprintf(stderr,"\x1B[0m");

#define DF2 \
t = clock() - t; \
fprintf(stderr,"\x1B[32m"); \
fprintf(stderr,"   ...%s ended (%f s lasted)\n\n",__PRETTY_FUNCTION__,((float)t)/CLOCKS_PER_SEC); \
fprintf(stderr,"\x1B[0m");
#else

#define DF1
#define DF2

#endif

#endif
