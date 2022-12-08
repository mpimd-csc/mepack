project(MEXOCT LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.0)
enable_testing()

#
# Options
#
option(EXAMPLES "Build Examples" ON)
option(MEXOCT_MATLAB_DOC     "Extract MATLAB documentation" ON)

include(${CMAKE_CURRENT_SOURCE_DIR}/MexOct.cmake)


if(EXAMPLES)
add_subdirectory(example)
endif()
add_subdirectory(tests)



