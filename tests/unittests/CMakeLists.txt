cmake_minimum_required(VERSION 2.6)
project(basic_test)

################################
# GTest
################################
add_subdirectory (gtest/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################
# Unit Tests
################################
# Add test cpp file
add_executable( runUnitTests test.cpp )
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main ulysseslib)
add_test( runUnitTests runUnitTests )
