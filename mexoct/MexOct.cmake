
SET(CMAKE_CXX_STANDARD 11)
INCLUDE(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++11" COMPILER_CXX_11)
if(NOT COMPILER_CXX_11)
    message(FATAL_ERROR "C++ compiler does not support C++11")
endif()

SET(MEXOCT_DIR ${CMAKE_CURRENT_LIST_DIR})
LIST(APPEND CMAKE_MODULE_PATH ${MEXOCT_DIR}/cmake)
INCLUDE_DIRECTORIES(${MEXOCT_DIR})

option(MEXOCT_MATLAB_DOC     "Extract MATLAB documentation" ON)
option(MEXOCT_MATLAB         "Search for Mathworks MATLAB" ON)
option(MEXOCT_OCTAVE         "Search for GNU Octave" ON)

if ( MEXOCT_MATLAB )
    find_package(MexOctMatlab)
endif()
if ( MEXOCT_OCTAVE )
    find_package(MexOctOctave)
endif()

if (NOT (MexOctOctave_FOUND OR MexOctMatlab_FOUND) )
    message(FATAL_ERROR "Octave and/or MATLAB missing. Building MEX/OCT interface not possible.")
endif()

if (MEXOCT_MATLAB_FOUND)
    mexoct_add_mex(extract_helptext ${MEXOCT_DIR}/tools/extract_helptext.c)
endif()

macro(mexoct_add target)
    if ( MEXOCT_MATLAB_FOUND )
        message(STATUS "Build MATLAB interface for ${target}")
        mexoct_add_mex(${target} ${ARGN})
        set_target_properties(mex_${target} PROPERTIES CXX_STANDARD 11)
        if (MEXOCT_MATLAB_DIRECTORY)
            set_target_properties(mex_${target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${MEXOCT_MATLAB_DIRECTORY}")
        else()
            set_target_properties(mex_${target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/matlab")

        endif()

        if (MEXOCT_MATLAB_DOC AND TARGET mex_extract_helptext)
            get_target_property(EXTRACT_HELP_NAME mex_extract_helptext OUTPUT_NAME)
            get_target_property(EXTRACT_HELP_DIR  mex_extract_helptext LIBRARY_OUTPUT_DIRECTORY)
            if ( NOT EXTRACT_HELP_DIR )
                get_target_property(EXTRACT_HELP_DIR  mex_extract_helptext BINARY_DIR)
            endif()

            get_target_property(C_HELP_NAME mex_${target} OUTPUT_NAME)
            get_target_property(C_HELP_DIR  mex_${target} LIBRARY_OUTPUT_DIRECTORY)
            if (NOT C_HELP_DIR)
                get_target_property(C_HELP_DIR  mex_${target} BINARY_DIR)
            endif()

            add_custom_target(mex_${target}_doc ALL COMMAND "${MEXOCT_MATLAB_ROOT}/bin/matlab" "-nojvm" "-nodesktop" "-nosplash"
                          -r "\"addpath('${EXTRACT_HELP_DIR}');extract_helptext('${C_HELP_DIR}/${C_HELP_NAME}${MEXOCT_MATLAB_EXT}','${C_HELP_DIR}/${C_HELP_NAME}.m');quit();\""
                          DEPENDS mex_${target} mex_extract_helptext
                          SOURCES ${ARGN}
                          BYPRODUCTS "${C_HELP_DIR}/${C_HELP_NAME}.m")
        endif()
    endif()

    if ( MEXOCT_OCTAVE_FOUND)
        message(STATUS "Build OCTAVE interface for ${target}")
        mexoct_add_oct(${target} ${ARGN})
        set_target_properties(oct_${target} PROPERTIES CXX_STANDARD 11)

        if (MEXOCT_OCTAVE_DIRECTORY)
            set_target_properties(oct_${target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${MEXOCT_OCTAVE_DIRECTORY}")
        else()
            set_target_properties(oct_${target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/octave")
        endif()
    endif()
endmacro()


macro(mexoct_matlab_target_link_libraries target)
    if ( MEXOCT_MATLAB_FOUND )
        target_link_libraries(mex_${target} ${ARGN})
    endif()
endmacro()

macro(mexoct_octave_target_link_libraries target)
    if ( MEXOCT_OCTAVE_FOUND)
        # foreach (LIB ${ARGN})
        #     message(STATUS ${LIB})
        #     if (TARGET ${LIB})
        #         MESSAGE("${LIB} is target")
        #         set_target_properties(${LIB} PROPERTIES LINK_LIBRARIES "" INTERFACE_LINK_LIBRARIES "")
        #     endif()
        # endforeach()
        target_link_libraries(oct_${target} ${ARGN})
    endif()
endmacro()

macro(mexoct_target_link_libraries target)
    mexoct_matlab_target_link_libraries(${target} ${ARGN})
    mexoct_octave_target_link_libraries(${target} ${ARGN})
endmacro()

macro(mexoct_target_link_blas target)
    mexoct_matlab_target_link_libraries(${target} ${MEXOCT_MATLAB_LAPACK_LIBRARY} ${MEXOCT_MATLAB_BLAS_LIBRARY})
    mexoct_octave_target_link_libraries(${target} ${MEXOCT_OCTAVE_LAPACK_LIBRARY} ${MEXOCT_OCTAVE_BLAS_LIBRARY})
endmacro()

macro(mexoct_target_include_directories target)
    if ( MEXOCT_MATLAB_FOUND )
        target_include_directories(mex_${target} ${ARGN})
    endif()
    if ( MEXOCT_OCTAVE_FOUND )
        target_include_directories(oct_${target} ${ARGN})
    endif()
endmacro()




macro(mexoct_install_target target)
    if ( MEXOCT_MATLAB_FOUND )
        install(TARGETS mex_${target} ${ARGN})
    endif()
    if ( MEXOCT_OCTAVE_FOUND )
        install(TARGETS oct_${target} ${ARGN})
    endif()
endmacro()
