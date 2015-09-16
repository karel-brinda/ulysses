cmake_minimum_required(VERSION 2.8)
#project(ulysses_tests)

################################
# GTest
################################

add_subdirectory (gtest/googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################
# Unit Tests
################################

add_executable(ulysses_tests ulysses_tests.cpp )
target_link_libraries(ulysses_tests gtest gtest_main ulysseslib ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
