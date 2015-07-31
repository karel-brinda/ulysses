cmake_minimum_required (VERSION 2.8)
project (ulysses)

#add_definitions(-O3 -W -Wall -Wextra -pedantic -std=c++11 -g -pg)
SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -W -Wall -Wextra -pedantic -std=c++11") # -g -pg")
SET(GCC_COVERAGE_LINK_FLAGS    "") #"-pg")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS       OFF) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost COMPONENTS regex serialization  REQUIRED)
include_directories("dynamic_bitset/include/"  ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})


add_executable(ulysses
    kmer.cpp
    kmer.hpp
    seed.cpp
    seed.hpp
    bloom.cpp
    bloom.hpp
    main.cpp
    main.hpp
    misc.hpp
    utils/MurmurHash3.cpp
    utils/MurmurHash3.hpp
    utils/kseq.h
)
target_link_libraries(ulysses m ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

