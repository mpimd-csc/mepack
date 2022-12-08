cmake_minimum_required(VERSION 3.5)
PROJECT(MEPACKTEST)
enable_language(Fortran)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(MEPACK)
find_package(fortran_stdlib REQUIRED)

ADD_EXECUTABLE(test test.c)
TARGET_LINK_LIBRARIES(test PRIVATE MEPACK::mepack)

ADD_EXECUTABLE(test_fortran test.f90)
TARGET_LINK_LIBRARIES(test_fortran PRIVATE MEPACK::mepack   fortran_stdlib::fortran_stdlib )


