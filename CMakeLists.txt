cmake_minimum_required (VERSION 2.8)

project (ulysses)
set(ULYSSES_MAJOR_VERSION 0)
set(ULYSSES_MINOR_VERSION 2)
set(ULYSSES_PATCH_VERSION 0)
set(ULYSSES_VERSION ${ULYSSES_MAJOR_VERSION}.${ULYSSES_MINOR_VERSION}.${ULYSSES_PATCH_VERSION})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

#add_definitions(-O3 -W -Wall -Wextra -pedantic -std=c++11 -g -pg)
SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -W -Wall -Wextra -pedantic -std=c++11") # -g -pg")
SET(GCC_COVERAGE_LINK_FLAGS    "") #"-pg")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )


include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")

include_directories("${PROJECT_SOURCE_DIR}/external/dynamic_bitset/include/")
include_directories("${PROJECT_SOURCE_DIR}/external/MurmurHash/")
include_directories("${PROJECT_SOURCE_DIR}/external/kseq")

#################

find_package(Boost COMPONENTS regex serialization  REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

set(Boost_USE_STATIC_LIBS       OFF) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)


#################

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIRS})

#################


# Add sub-directories
add_subdirectory(lib)
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")


add_executable(ulysses
    main/main.cpp
    main/main.hpp
    )

target_link_libraries(ulysses
        m
        ${ZLIB_LIBRARIES}
        ${Boost_LIBRARIES}
        ulysseslib
    )
