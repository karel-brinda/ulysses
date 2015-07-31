//
//  seed.c
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//
//

#include "seed.hpp"

seed_t::seed_t(const char *seedstr) {
    seed_t * seed = this;
    
    if(seedstr==NULL){
        seed->span=0;
        seed->weight=0;
        return ;
    }
    
    for(int i=0;i<MAX_SEED_SPAN;i++){
        seed->seedstr[i]=0xa;
        seed->shape[0][i]=0xa;
        seed->shape[1][i]=0xa;
    }
        
    int span=strlen(seedstr);
    assert(span<MAX_SEED_SPAN);
    strcpy(seed->seedstr,seedstr);
    
    int w=0;
    for(int i=0;i<span;i++){
        assert(seedstr[i]=='#' || seedstr[i]=='-');
        if (seedstr[i]=='#'){
            w++;
        }
    }

    seed->span=span;
    seed->weight=w;

    int v=0;
    for(int i=0;i<span;i++)
    {
        if (seedstr[i]=='#'){
            seed->shape[0][v]=i;
            seed->shape[1][v]=span-i-1;
            v++;
        }
    }
    
    return ;
}
