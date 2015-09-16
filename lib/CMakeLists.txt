cmake_minimum_required(VERSION 2.8)

add_library (ulysseslib 
        kmer.cpp
        seed.cpp
        bloom.cpp
        ../external/MurmurHash/MurmurHash3.cpp
    )
