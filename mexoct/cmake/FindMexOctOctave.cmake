# - Find Octave
# GNU Octave is a high-level interpreted language, primarily intended for numerical computations.
# available at http://www.gnu.org/software/octave/
#
# This module defines:
#  OCTAVE_EXECUTABLE           - octave interpreter
#  OCTAVE_MKOCTFILE            - mkoctfile command
#  OCTAVE_INCLUDE_DIRS         - include path for mex.h, mexproto.h
#  OCTAVE_LIBRARIES            - required libraries: octinterp, octave, cruft
#  OCTAVE_OCTINTERP_LIBRARY    - path to the library octinterp
#  OCTAVE_OCTAVE_LIBRARY       - path to the library octave
#  OCTAVE_CRUFT_LIBRARY        - path to the library cruft
#  OCTAVE_VERSION_STRING       - octave version string
#  OCTAVE_MAJOR_VERSION        - major version
#  OCTAVE_MINOR_VERSION        - minor version
#  OCTAVE_PATCH_VERSION        - patch version
#  OCTAVE_OCT_FILE_DIR         - object files that will be dynamically loaded
#  OCTAVE_OCT_LIB_DIR          - oct libraries
#  OCTAVE_ROOT_DIR             - octave prefix
#
# The macro octave_add_oct allows to create compiled modules.
# octave_add_oct ( target_name function_name
#         [SOURCES] source1 [source2 ...]
#         [LINK_LIBRARIES  lib1 [lib2 ...]]
#         [EXTENSION ext]
# )
#
# To install it, you can the use the variable OCTAVE_OCT_FILE_DIR as follow:
#  file ( RELATIVE_PATH PKG_OCTAVE_OCT_FILE_DIR ${OCTAVE_ROOT_DIR} ${OCTAVE_OCT_FILE_DIR} )
#  install (
#    TARGETS target_name
#    DESTINATION ${PKG_OCTAVE_OCT_FILE_DIR}
#  )

#=============================================================================
# Copyright 2013, Julien Schueller
# Copyright 2015, Martin Koehler
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of the FreeBSD Project.
#=============================================================================

if ( NOT MEXOCT_OCTAVE_CONFIG )
    find_program(MEXOCT_OCTAVE_CONFIG_EXECUTABLE
        NAMES octave-config
    )
else()
    set (MEXOCT_OCTAVE_CONFIG_EXECUTABLE ${MEXOCT_OCTAVE_CONFIG})
endif()

MESSAGE(STATUS "Octave config comes from ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE}")
MESSAGE(STATUS "A user defined one can be supplied by MEXOCT_OCTAVE_CONFIG")

if( MEXOCT_OCTAVE_CONFIG_EXECUTABLE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p PREFIX
        OUTPUT_VARIABLE MEXOCT_OCTAVE_ROOT_DIR1
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p OCTAVE_HOME
        OUTPUT_VARIABLE MEXOCT_OCTAVE_ROOT_DIR2
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p BINDIR
        OUTPUT_VARIABLE MEXOCT_OCTAVE_BIN_PATHS
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p OCTINCLUDEDIR
        OUTPUT_VARIABLE MEXOCT_OCTAVE_INCLUDE_PATHS
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p OCTLIBDIR
        OUTPUT_VARIABLE MEXOCT_OCTAVE_LIBRARIES_PATHS
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p OCTFILEDIR
        OUTPUT_VARIABLE MEXOCT_OCTAVE_OCT_FILE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -p OCTLIBDIR
        OUTPUT_VARIABLE MEXOCT_OCTAVE_OCT_LIB_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    execute_process( COMMAND ${MEXOCT_OCTAVE_CONFIG_EXECUTABLE} -v
        OUTPUT_VARIABLE MEXOCT_OCTAVE_VERSION_STRING
        OUTPUT_STRIP_TRAILING_WHITESPACE )

    if( MEXOCT_OCTAVE_ROOT_DIR1 )
        set (MEXOCT_OCTAVE_ROOT_DIR ${MEXOCT_OCTAVE_ROOT_DIR1})
    else()
        set (MEXOCT_OCTAVE_ROOT_DIR ${MEXOCT_OCTAVE_ROOT_DIR2})
    endif()

    if( MEXOCT_OCTAVE_VERSION_STRING )
        string( REGEX REPLACE "([0-9]+)\\..*" "\\1" MEXOCT_OCTAVE_MAJOR_VERSION ${MEXOCT_OCTAVE_VERSION_STRING} )
        string( REGEX REPLACE "[0-9]+\\.([0-9]+).*" "\\1" MEXOCT_OCTAVE_MINOR_VERSION ${MEXOCT_OCTAVE_VERSION_STRING} )
        string( REGEX REPLACE "[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" MEXOCT_OCTAVE_PATCH_VERSION ${MEXOCT_OCTAVE_VERSION_STRING} )
    endif()
endif()

find_program(MEXOCT_OCTAVE_EXECUTABLE
    NAMES octave octave-${MEXOCT_OCTAVE_MAJOR_VERSION}.${MEXOCT_OCTAVE_MINOR_VERSION}.${MEXOCT_OCTAVE_PATCH_VERSION}
    HINTS ${MEXOCT_OCTAVE_BIN_PATHS}
    )

find_program(MEXOCT_OCTAVE_MKOCTFILE
    NAMES mkoctfile mkoctfile-${MEXOCT_OCTAVE_MAJOR_VERSION}.${MEXOCT_OCTAVE_MINOR_VERSION}.${MEXOCT_OCTAVE_PATCH_VERSION} mkoctfile
    HINTS ${MEXOCT_OCTAVE_BIN_PATHS}
    )


find_library(MEXOCT_OCTAVE_OCTINTERP_LIBRARY
    NAMES octinterp liboctinterp
    HINTS ${MEXOCT_OCTAVE_LIBRARIES_PATHS}
    )
find_library(MEXOCT_OCTAVE_LIBRARY
    NAMES octave liboctave
    HINTS ${MEXOCT_OCTAVE_LIBRARIES_PATHS}
    )
find_library(MEXOCT_OCTAVE_CRUFT_LIBRARY
    NAMES cruft libcruft
    HINTS ${MEXOCT_OCTAVE_LIBRARIES_PATHS}
    )

set( MEXOCT_OCTAVE_LIBRARIES ${MEXOCT_OCTAVE_OCTINTERP_LIBRARY} )
list( APPEND MEXOCT_OCTAVE_LIBRARIES ${MEXOCT_OCTAVE_LIBRARY} )
if(${MEXOCT_OCTAVE_CRUFT_LIBRARY})
    list( APPEND MEXOCT_OCTAVE_LIBRARIES ${MEXOCT_OCTAVE_CRUFT_LIBRARY} )
endif()

find_path(MEXOCT_OCTAVE_INCLUDE_DIR
    NAMES oct.h
    HINTS ${MEXOCT_OCTAVE_INCLUDE_PATHS}
    )
execute_process(COMMAND dirname ${MEXOCT_OCTAVE_INCLUDE_DIR} OUTPUT_VARIABLE MEXOCT_OCTAVE_INCLUDE_BASE OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process( COMMAND ${MEXOCT_OCTAVE_MKOCTFILE} -p BLAS_LIBS
        OUTPUT_VARIABLE MEXOCT_OCTAVE_BLAS_LIBRARY
        OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND ${MEXOCT_OCTAVE_MKOCTFILE} -p LAPACK_LIBS
        OUTPUT_VARIABLE MEXOCT_OCTAVE_LAPACK_LIBRARY
        OUTPUT_STRIP_TRAILING_WHITESPACE )

set( MEXOCT_OCTAVE_INCLUDE_DIR ${MEXOCT_OCTAVE_INCLUDE_DIR} ${MEXOCT_OCTAVE_INCLUDE_BASE} )

macro( mexoct_add_oct target )
    set( _CMD SOURCES )
    set( _SOURCES )
    set( _LINK_LIBRARIES )
    set( _EXT_INCLUDE_DIRS )
    set( _EXTENSION )
    set( _OCT_EXTENSION oct )
    foreach(_ARG ${ARGN})
        if( ${_ARG} MATCHES SOURCES )
            set( _CMD SOURCES )
        elseif( ${_ARG} MATCHES LINK_LIBRARIES )
            set( _CMD LINK_LIBRARIES )
        elseif( ${_ARG} MATCHES EXT_INCLUDE_DIRS )
            set( _CMD EXT_INCLUDE_DIRS )
        elseif( ${_ARG} MATCHES EXTENSION )
            set( _CMD EXTENSION )
        else()
            if( ${_CMD} MATCHES SOURCES )
                list( APPEND _SOURCES "${_ARG}" )
            elseif(${_CMD} MATCHES LINK_LIBRARIES)
                list( APPEND _LINK_LIBRARIES "${_ARG}" )
            elseif(${_CMD} MATCHES EXT_INCLUDE_DIRS)
                list( APPEND _EXT_INCLUDE_DIRS "${_ARG}" )
            elseif( ${_CMD} MATCHES EXTENSION )
                set( _OCT_EXTENSION ${_ARG} )
            endif()
        endif()
    endforeach()

    if(NOT ${target}_NAME)
        set(${target}_NAME oct_${target})
    endif()

    if(NOT ${target}_OUTPUT_NAME)
        set(${target}_OUTPUT_NAME ${target})
    endif()

    # include_directories(${MEXOCT_OCTAVE_INCLUDE_DIRS} ${_EXT_INCLUDE_DIRS})
    add_library( ${${target}_NAME} SHARED ${_SOURCES} )
    target_link_libraries( ${${target}_NAME} ${MEXOCT_OCTAVE_LIBRARIES} ${_LINK_LIBRARIES} )
    target_include_directories( ${${target}_NAME} PRIVATE ${MEXOCT_OCTAVE_INCLUDE_DIR} ${_EXT_INCLUDE_DIRS})
    target_compile_definitions( ${${target}_NAME} PRIVATE "MEXOCT_OCTAVE")
    set_target_properties(${${target}_NAME} PROPERTIES
        PREFIX ""
        SUFFIX  ".${_OCT_EXTENSION}"
        OUTPUT_NAME "${target}"
        )
endmacro()

# handle REQUIRED and QUIET options
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args(MexOctOctave REQUIRED_VARS MEXOCT_OCTAVE_EXECUTABLE MEXOCT_OCTAVE_ROOT_DIR MEXOCT_OCTAVE_INCLUDE_DIR MEXOCT_OCTAVE_INCLUDE_PATHS MEXOCT_OCTAVE_LIBRARIES VERSION_VAR MEXOCT_OCTAVE_VERSION_STRING)

IF ( MexOctOctave_FOUND )
    SET(MEXOCT_OCTAVE_FOUND TRUE)
ENDIF()

mark_as_advanced(
    MEXOCT_OCTAVE_OCT_FILE_DIR
    MEXOCT_OCTAVE_OCT_LIB_DIR
    MEXOCT_OCTAVE_OCTINTERP_LIBRARY
    MEXOCT_OCTAVE_LIBRARY
    MEXOCT_OCTAVE_CRUFT_LIBRARY
    MEXOCT_OCTAVE_LIBRARIES
    MEXOCT_OCTAVE_INCLUDE_DIR
    MEXOCT_OCTAVE_INCLUDE_DIRS
    MEXOCT_OCTAVE_INCLUDE_BASE
    MEXOCT_OCTAVE_ROOT_DIR
    MEXOCT_OCTAVE_VERSION_STRING
    MEXOCT_OCTAVE_MAJOR_VERSION
    MEXOCT_OCTAVE_MINOR_VERSION
    MEXOCT_OCTAVE_PATCH_VERSION
    )
