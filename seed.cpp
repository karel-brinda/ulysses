//
//  seed.c
//  bloom
//
//  Created by Karel Brinda on 29.04.15.
//
//

#include "seed.hpp"

Seed::Seed(const Seed &seed)
    :weight(seed.weight),span(seed.span)
{
    strcpy(this->seedstr,seed.seedstr);

    memcpy(&this->shape[0], seed.shape[0], MAX_SEED_SPAN);
    memcpy(&this->shape[1], seed.shape[1], MAX_SEED_SPAN);
}


Seed::Seed(const char *seedstr) {
    
    if(seedstr==NULL){
        this->span=0;
        this->weight=0;
        return ;
    }
    
    for(int i=0;i<MAX_SEED_SPAN;i++){
        this->seedstr[i]=0xa;
        this->shape[0][i]=0xa;
        this->shape[1][i]=0xa;
    }
        
    int span=strlen(seedstr);
    assert(span<MAX_SEED_SPAN);
    strcpy(this->seedstr,seedstr);
    
    int w=0;
    for(int i=0;i<span;i++){
        assert(seedstr[i]=='#' || seedstr[i]=='-');
        if (seedstr[i]=='#'){
            w++;
        }
    }

    this->span=span;
    this->weight=w;

    int v=0;
    for(int i=0;i<span;i++)
    {
        if (seedstr[i]=='#'){
            this->shape[0][v]=i;
            this->shape[1][v]=span-i-1;
            v++;
        }
    }
    
    return ;
}
