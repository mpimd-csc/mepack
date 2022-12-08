PROJECT(MEXOCT_DEMO)
cmake_minimum_required(VERSION 3.0)

INCLUDE(${PROJECT_SOURCE_DIR}/../MexOct.cmake)

IF (MEXOCT_MATLAB_FOUND)
    mexoct_add(lu lu.cpp)
ENDIF()

