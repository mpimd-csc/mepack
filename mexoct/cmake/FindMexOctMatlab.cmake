# - this module looks for Matlab
# Defines:
#  MEXOCT_MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MEXOCT_MATLAB_LIBRARIES:   required libraries: libmex, etc
#  MEXOCT_MATLAB_MEX_LIBRARY: path to libmex.lib
#  MEXOCT_MATLAB_MX_LIBRARY:  path to libmx.lib
#  MEXOCT_MATLAB_ENG_LIBRARY: path to libeng.lib
#  MEXOCT_MATLAB_ROOT: path to Matlab's root directory
#  MEXOCT_MATLAB_ARCH_DIR:    MATLAB ARCH DIR
#  MEXOCT_MATLAB_EXT


# Get all propreties that cmake supports
execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

# Convert command output into a CMake list
STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

function(print_properties)
    message ("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
endfunction(print_properties)

function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
        message("There is no target named '${tgt}'")
        return()
    endif()

    foreach (prop ${CMAKE_PROPERTY_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
        # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
        if(prop STREQUAL "LOCATION" OR prop MATCHES "^LOCATION_" OR prop MATCHES "_LOCATION$")
            continue()
        endif()
        # message ("Checking ${prop}")
        get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
        if (propval)
            get_target_property(propval ${tgt} ${prop})
            message ("${tgt} ${prop} = ${propval}")
        endif()
    endforeach(prop)
endfunction(print_target_properties)

function(get_value_from_xml INPUT TAG REG OUTPUT)
    STRING(FIND "${INPUT}" "<${TAG}>" START)
    STRING(FIND "${INPUT}" "</${TAG}>" STOP)
    MATH(EXPR LENGTH "${STOP} - ${START}")
    STRING(SUBSTRING "${INPUT}" "${START}" "${LENGTH}" DATA)
    IF(DATA MATCHES "^<${TAG}>${REG}$")
        SET(${OUTPUT} ${CMAKE_MATCH_1} PARENT_SCOPE)
    ELSE()
        SET(${OUTPUT} "NOTFOUND" PARENT_SCOPE)
    ENDIF()
endfunction()

INCLUDE(FindPackageHandleStandardArgs)
INCLUDE(CheckCXXCompilerFlag)

SET(MEXOCT_MATLAB_FOUND 0)
IF(WIN32)
    # Search for a version of Matlab available, starting from the most modern one to older versions
    FOREACH(MATVER  "24.2" "24.1"
            "23.2"
            "9.14" "9.13" "9.12" "9.11" "9.10" "9.9" "9.8" "9.7" "9.6" "9.5" "9.4" "9.3" "9.2" "9.1" "9.0"
            "8.6" "8.5" "8.4" "8.3" "8.2" "8.1"
            "7.14" "7.13" "7.12" "7.11" "7.10" "7.9" "7.8" "7.7" "7.6" "7.5" "7.4")
        IF((NOT DEFINED MEXOCT_MATLAB_ROOT) OR ("${MEXOCT_MATLAB_ROOT}" STREQUAL "") OR("${MEXOCT_MATLAB_ROOT}" STREQUAL "/registry"))
            GET_FILENAME_COMPONENT(MEXOCT_MATLAB_ROOT
                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\${MATVER};MATLABROOT]"
                ABSOLUTE)
            SET(MEXOCT_MATLAB_VERSION ${MATVER})
        ENDIF()
    ENDFOREACH()

    IF(NOT EXISTS "${MEXOCT_MATLAB_ROOT}")
        MESSAGE(WARNING "MEXOCT_MATLAB_ROOT = ${MEXOCT_MATLAB_ROOT} does not exist. MATLAB could not be used.")
        RETURN()
    ENDIF()

    # Directory name depending on whether the Windows architecture is 32
    # bit or 64 bit
    IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
        SET(WINDIR "win32")
        SET(MEXOCT_MATLAB_EXT ".mexw32")
    ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET(WINDIR "win64")
        SET(MEXOCT_MATLAB_EXT ".mexw64")
    ELSE()
        MESSAGE(FATAL_ERROR "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
    ENDIF()

    # Folder where the MEX libraries are, depending of the Windows compiler
    IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
        SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/msvc60")
    ELSEIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
        # Assume people are generally using Visual Studio 7.1,
        # if using 7.0 need to link to: ../extern/lib/${WINDIR}/microsoft/msvc70
        SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/msvc71")
        # SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/msvc70")
    ELSEIF(${CMAKE_GENERATOR} MATCHES "Borland")
        # Assume people are generally using Borland 5.4,
        # if using 7.0 need to link to: ../extern/lib/${WINDIR}/microsoft/msvc70
        SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/bcc54")
        # SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/bcc50")
        # SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/bcc51")
    ELSEIF(${CMAKE_GENERATOR} MATCHES "Visual Studio*")
        # If the compiler is Visual Studio, but not any of the specific
        # versions above, we try our luck with the microsoft directory
        SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/microsoft/")
    ELSEIF(MINGW AND CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET(MEXOCT_MATLAB_LIBRARIES_DIR "${MEXOCT_MATLAB_ROOT}/extern/lib/${WINDIR}/mingw64")
    ELSE()
        MESSAGE(FATAL_ERROR "Generator not compatible: ${CMAKE_GENERATOR}")
    ENDIF()

    MESSAGE(STATUS "MATLAB ROOT = ${MEXOCT_MATLAB_ROOT}")

    IF (EXISTS "${MEXOCT_MATLAB_ROOT}/VersionInfo.xml")
        FILE(READ "${MEXOCT_MATLAB_ROOT}/VersionInfo.xml" VERSION_INFO_FILE)
        get_value_from_xml("${VERSION_INFO_FILE}" "version" "([0-9]+\.[0-9]+)\..*" MEXOCT_MATLAB_VERSION)
        get_value_from_xml("${VERSION_INFO_FILE}" "release" "([0-9a-zA-Z]+)" MEXOCT_MATLAB_RELEASE)
    ENDIF()


    # Get paths to the Matlab MEX libraries
    FIND_LIBRARY(MEXOCT_MATLAB_MEX_LIBRARY
        libmex
        ${MEXOCT_MATLAB_LIBRARIES_DIR}
        )
    FIND_LIBRARY(MEXOCT_MATLAB_MX_LIBRARY
        libmx
        ${MEXOCT_MATLAB_LIBRARIES_DIR}
        )
    FIND_LIBRARY(MEXOCT_MATLAB_ENG_LIBRARY
        libeng
        ${MEXOCT_MATLAB_LIBRARIES_DIR}
        )
    FIND_LIBRARY(MEXOCT_MATLAB_BLAS_LIBRARY
        libmwblas
        ${MEXOCT_MATLAB_LIBRARIES_DIR}
        )
    FIND_LIBRARY(MEXOCT_MATLAB_LAPACK_LIBRARY
        libmwlapack
        ${MEXOCT_MATLAB_LIBRARIES_DIR}
        )
    # Get path to the include directory
    FIND_PATH(MEXOCT_MATLAB_INCLUDE_DIR
        "mex.h"
        "${MEXOCT_MATLAB_ROOT}\\extern\\include"
        )
    FIND_FILE(MEXOCT_MATLAB_MEX_VERSION_FILE
        NAMES c_mexapi_version.c
        PATHS "${MEXOCT_MATLAB_ROOT}\\extern\\version")
    MESSAGE(STATUS "MEXOCT_MATLAB_MEX_VERSION_FILE = ${MEXOCT_MATLAB_MEX_VERSION_FILE}")

    FIND_PROGRAM(MEXOCT_MATLAB_BINARY NAMES matlab.exe PATHS "${MEXOCT_MATLAB_ROOT}\\bin" NO_DEFAULT_PATH)
    FIND_PROGRAM(MEXOCT_MEX_BINARY NAMES mex mex.exe mex.bat PATHS "${MEXOCT_MATLAB_ROOT}\\bin" NO_DEFAULT_PATH)

    SET(MEXOCT_MATLAB_BINARY ${MEXOCT_MATLAB_BINARY} CACHE FILEPATH "MATLAB Executable")
    SET(MEXOCT_MEX_BINARY ${MEXOCT_MEX_BINARY} CACHE FILEPATH "MATLAB MEX Executable")


    MESSAGE(STATUS "Found MATLAB (${MEXOCT_MATLAB_RELEASE} - ${MEXOCT_MATLAB_VERSION})")

ELSE()
    IF ( DEFINED MEXOCT_MATLAB_ROOT AND NOT "${MEXOCT_MATLAB_ROOT}" STREQUAL "")


    ELSE ()
        IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
            # Search for a version of Matlab available, starting from the most modern one to older versions
            FOREACH(MATVER "R2025b" "R2025a" "R2024b" "R2024a" "R2023b" "R2023a" "R2022b" "R2022a" "R2021b" "R2021a" "R2020b" "R2020a"
                    "R2019b" "R2019a" "R2018b" "R2018a" "R2017b" "R2017a" "R2016b" "R2016a" "R2015b" "R2015a" "R2014b" "R2014a"
                    "R2013b" "R2013a" "R2012b" "R2012a" "R2011b" "R2011a" "R2010b" "R2010a" "R2009b" "R2009a" "R2008b")
                IF((NOT DEFINED MEXOCT_MATLAB_ROOT) OR ("${MEXOCT_MATLAB_ROOT}" STREQUAL ""))
                    IF(EXISTS /Applications/MEXOCT_MATLAB_${MATVER}.app)
                        SET(MATLAB_ROOT /Applications/MATLAB_${MATVER}.app)

                    ENDIF()
                ENDIF()
            ENDFOREACH()

        ELSE()
            EXECUTE_PROCESS(
                COMMAND which matlab
                COMMAND xargs readlink
                COMMAND xargs dirname
                COMMAND xargs dirname
                COMMAND xargs echo -n
                OUTPUT_VARIABLE MEXOCT_MATLAB_ROOT
                )
        ENDIF()
    ENDIF()
    SET(MEXOCT_MATLAB_ROOT ${MEXOCT_MATLAB_ROOT} CACHE PATH "MATLAB Root Directory")

    FIND_PROGRAM(MEXOCT_MATLAB_BINARY NAMES matlab PATHS ${MEXOCT_MATLAB_ROOT}/bin NO_DEFAULT_PATH)
    FIND_PROGRAM(MEXOCT_MEX_BINARY NAMES mex PATHS ${MEXOCT_MATLAB_ROOT}/bin NO_DEFAULT_PATH)

    SET(MEXOCT_MATLAB_BINARY ${MEXOCT_MATLAB_BINARY} CACHE FILEPATH "MATLAB Executable")
    SET(MEXOCT_MEX_BINARY ${MEXOCT_MEX_BINARY} CACHE FILEPATH "MATLAB MEX Executable")


    IF ( MEXOCT_MATLAB_BINARY AND MEXOCT_MEX_BINARY )

        MESSAGE(STATUS "MEXOCT_MATLAB_ROOT = ${MEXOCT_MATLAB_ROOT}")
        # Check if this is a Mac
        IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
            SET(LIBRARY_EXTENSION .dylib)
        ELSE()
            SET(LIBRARY_EXTENSION .so)
        ENDIF()

        # Get path to the MEX libraries
        IF(APPLE)
            IF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
                SET(MEXOCT_MATLAB_EXT ".mexmaci")
            ELSE()
                SET(MEXOCT_MATLAB_EXT ".mexmac")
            ENDIF()
        ELSE()
            IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
                SET(MEXOCT_MATLAB_EXT ".mexglx")
                SET(MEXOCT_MATLAB_ARCH_DIR "glnx86")
            ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
                SET(MEXOCT_MATLAB_EXT ".mexa64")
                SET(MEXOCT_MATLAB_ARCH_DIR "glnxa64")
            ELSE()
                MESSAGE(FATAL_ERROR "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
            ENDIF()
        ENDIF()

        MESSAGE(STATUS "MEXOCT_MATLAB_EXT = ${MEXOCT_MATLAB_EXT}")
        SET(MEXOCT_MATLAB_EXT ${MEXOCT_MATLAB_EXT} CACHE STRING "MATLAB Mex Extension")

        SET(MEXOCT_MATLAB_ARCH_DIR "${MEXOCT_MATLAB_ROOT}/bin/${MEXOCT_MATLAB_ARCH_DIR}/" CACHE PATH "MATLAB Arch Directory")

        IF(EXISTS ${MEXOCT_MATLAB_ARCH_DIR})
            EXECUTE_PROCESS(
                COMMAND find "${MEXOCT_MATLAB_ARCH_DIR}" -name libmex${LIBRARY_EXTENSION}
                COMMAND xargs echo -n
                OUTPUT_VARIABLE MEXOCT_MATLAB_MEX_LIBRARY
                )
            EXECUTE_PROCESS(
                COMMAND find "${MEXOCT_MATLAB_ARCH_DIR}" -name libmx${LIBRARY_EXTENSION}
                COMMAND xargs echo -n
                OUTPUT_VARIABLE MEXOCT_MATLAB_MX_LIBRARY
                )
            EXECUTE_PROCESS(
                COMMAND find "${MEXOCT_MATLAB_ARCH_DIR}" -name libeng${LIBRARY_EXTENSION}
                COMMAND xargs echo -n
                OUTPUT_VARIABLE MEXOCT_MATLAB_ENG_LIBRARY
                )
            EXECUTE_PROCESS(
                COMMAND find "${MEXOCT_MATLAB_ARCH_DIR}" -name libmwblas${LIBRARY_EXTENSION}
                COMMAND xargs echo -n
                OUTPUT_VARIABLE MEXOCT_MATLAB_BLAS_LIBRARY
                )
            EXECUTE_PROCESS(
                COMMAND find "${MEXOCT_MATLAB_ARCH_DIR}" -name libmwlapack${LIBRARY_EXTENSION}
                COMMAND xargs echo -n
                OUTPUT_VARIABLE MEXOCT_MATLAB_LAPACK_LIBRARY
                )
            SET(MEXOCT_MATLAB_MEX_LIBRARY ${MEXOCT_MATLAB_MEX_LIBRARY} CACHE INTERNAL "Mex Library")
            SET(MEXOCT_MATLAB_MX_LIBRARY ${MEXOCT_MATLAB_MX_LIBRARY} CACHE INTERNAL "MX Library")
            SET(MEXOCT_MATLAB_ENG_LIBRARY ${MEXOCT_MATLAB_ENG_LIBRARY} CACHE INTERNAL "ENG Library")
            SET(MEXOCT_MATLAB_BLAS_LIBRARY ${MEXOCT_MATLAB_BLAS_LIBRARY} CACHE INTERNAL "mwblas Library")
            SET(MEXOCT_MATLAB_LAPACK_LIBRARY ${MEXOCT_MATLAB_LAPACK_LIBRARY} CACHE INTERNAL "mwlapack Library")
        ENDIF()



        # Get path to the include directory
        FIND_FILE(MEXOCT_MATLAB_MEX_VERSION_FILE
            NAMES c_mexapi_version.c
            PATHS "${MEXOCT_MATLAB_ROOT}/extern/version")
        MESSAGE(STATUS "MEXOCT_MATLAB_MEX_VERSION_FILE = ${MEXOCT_MATLAB_MEX_VERSION_FILE}")

        FIND_PATH(MEXOCT_MATLAB_INCLUDE_DIR
            "mex.h"
            PATHS "${MEXOCT_MATLAB_ROOT}/extern/include"
            )
        SET(MEXOCT_MATLAB_MEX ${MEXOCT_MATLAB_ROOT}/bin/mex CACHE FILEPATH "Final MEX executable")
        MESSAGE(STATUS "MEXOCT_MATLAB_MEX = ${MEXOCT_MATLAB_MEX}")



        # determine if matlab version is 9.5 or higher
        # https://de.mathworks.com/help/matlab/matlab_external/do-i-need-to-upgrade-my-mex-files-to-use-interleaved-complex.html
        IF (NOT DEFINED MEXOCT_MATLAB_MEX_INTERLEAVED_COMPLEX_API)

            EXECUTE_PROCESS(
                COMMAND "${MEXOCT_MATLAB_ROOT}/bin/matlab" -nosplash -nojvm -nodesktop -nodisplay -r "exit(~verLessThan('matlab','9.4'))"
                RESULT_VARIABLE  MEXOCT_MATLAB_MEX_INTERLEAVED_COMPLEX_API
                ERROR_VARIABLE  __OUTPUTE
                OUTPUT_VARIABLE __OUTPUTS
                TIMEOUT 180
                )
            IF (MEXOCT_MATLAB_MEX_INTERLEAVED_COMPLEX_API)
                MESSAGE( STATUS "MEXOCT_MATLAB can use the INTERLEAVED COMPLEX API")
            ENDIF()
            SET(MEXOCT_MATLAB_MEX_INTERLEAVED_COMPLEX_API ${MEXOCT_MATLAB_MEX_INTERLEAVED_COMPLEX_API} CACHE INTERNAL "Interleaved MEX Complex API")
        ENDIF()

        IF (NOT DEFINED MEXOCT_MATLAB_VERSION)

            EXECUTE_PROCESS(
                COMMAND "${MEXOCT_MATLAB_ROOT}/bin/matlab" -nodesktop -nojvm -nosplash -r "f=fopen('${CMAKE_CURRENT_BINARY_DIR}/matlab_version.txt','w');x=ver('matlab'); fprintf(f,'%s %s',x.Version,x.Release); fclose(f); quit()"
                RESULT_VARIABLE _X
                ERROR_VARIABLE  __OUTPUTE
                OUTPUT_VARIABLE __OUTPUTS
                TIMEOUT 180
                )


            FILE(READ "${CMAKE_CURRENT_BINARY_DIR}/matlab_version.txt" MEXOCT_MATLAB_VERSION_STRING)
            STRING(REPLACE " " ";" MEXOCT_MATLAB_VERSION_LIST ${MEXOCT_MATLAB_VERSION_STRING})
            LIST(GET MEXOCT_MATLAB_VERSION_LIST 0 MEXOCT_MATLAB_VERSION)
            LIST(GET MEXOCT_MATLAB_VERSION_LIST 1 MEXOCT_MATLAB_RELEASE)
            STRING(REPLACE "(" "" MEXOCT_MATLAB_RELEASE ${MEXOCT_MATLAB_RELEASE})
            STRING(REPLACE ")" "" MEXOCT_MATLAB_RELEASE ${MEXOCT_MATLAB_RELEASE})
            SET(MEXOCT_MATLAB_VERSION ${MEXOCT_MATLAB_VERSION} CACHE INTERNAL "MATLAB Version")
            SET(MEXOCT_MATLAB_RELEASE ${MEXOCT_MATLAB_RELEASE} CACHE INTERNAL "MATLAB Release")
        ENDIF()


        MESSAGE(STATUS "MEXOCT_MATLAB_VERSION is ${MEXOCT_MATLAB_VERSION}")
        MESSAGE(STATUS "MEXOCT_MATLAB_RELEASE is ${MEXOCT_MATLAB_RELEASE}")

    ELSE()
        MESSAGE(STATUS "Matlab or Mex Program not found")
    ENDIF()

ENDIF()

# This is common to UNIX and Win32:
SET(MEXOCT_MATLAB_LIBRARIES
    ${MEXOCT_MATLAB_MEX_LIBRARY}
    ${MEXOCT_MATLAB_MX_LIBRARY}
    )


IF(MEXOCT_MATLAB_INCLUDE_DIR AND MEXOCT_MATLAB_LIBRARIES AND MEXOCT_MATLAB_BINARY AND MEX_BINARY)
    SET(MEXOCT_MATLAB_FOUND 1)
ENDIF()

MARK_AS_ADVANCED(
    MEXOCT_MATLAB_LIBRARIES
    MEXOCT_MATLAB_MEX_LIBRARY
    MEXOCT_MATLAB_MX_LIBRARY
    MEXOCT_MATLAB_ENG_LIBRARY
    MEXOCT_MATLAB_INCLUDE_DIR
    MEXOCT_MATLAB_FOUND
    MEXOCT_MATLAB_ROOT
    )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MexOctMatlab  DEFAULT_MSG  MEXOCT_MATLAB_LIBRARIES MEXOCT_MATLAB_INCLUDE_DIR MEXOCT_MEX_BINARY MEXOCT_MATLAB_BINARY)
IF ( MexOctMatlab_FOUND )
    SET ( MEXOCT_MATLAB_FOUND TRUE)
ENDIF()

function(mexoct_add_mex target)
    if(NOT WIN32)
        # we do not need all this on Windows
        # pthread options
        check_cxx_compiler_flag(-pthread                    HAS_MINUS_PTHREAD)
        # we should use try_compile instead, the link flags are discarded from
        # this compiler_flag function.
        #check_cxx_compiler_flag(-Wl,--exclude-libs,ALL HAS_SYMBOL_HIDING_CAPABILITY)

    endif()

    # set(oneValueArgs NAME DOCUMENTATION OUTPUT_NAME)
    # set(multiValueArgs LINK_TO SRC)
    #
    # set(prefix _matlab_addmex_prefix)
    # cmake_parse_arguments(${prefix} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
    #

    if(NOT ${target}_NAME)
        set(${target}_NAME mex_${target})
    endif()

    if(NOT ${target}_OUTPUT_NAME)
        set(${target}_OUTPUT_NAME ${target})
    endif()

    if(MEXOCT_MATLAB_MEX_VERSION_FILE)
        add_library(${${target}_NAME} SHARED ${ARGN} ${MEXOCT_MATLAB_MEX_VERSION_FILE})
    else()
        add_library(${${target}_NAME} SHARED ${ARGN})
    endif()

    target_include_directories(${${target}_NAME} PRIVATE ${MEXOCT_MATLAB_INCLUDE_DIR})

    if(DEFINED MEXOCT_MATLAB_MX_LIBRARY)
        target_link_libraries(${${target}_NAME} ${MEXOCT_MATLAB_MX_LIBRARY})
    endif()

    if (MATLAB_MEX_INTERLEAVED_COMPLEX_API)
        set(MEX_API_MACRO "MATLAB_DEFAULT_RELEASE=R2018a")
        set(MEX_COMPLEX_INTERLEAVE "MEX_COMPLEX_INTERLEAVE=1")
        message(STATUS "Build mexfile ${target} with R2018a API")
    else()
        set(MEX_API_MACRO "MATLAB_DEFAULT_RELEASE=R2017b")
        set(MEX_COMPLEX_INTERLEAVE "MEX_COMPLEX_INTERLEAVE=0")
    endif()

    target_link_libraries(${${target}_NAME} ${MEXOCT_MATLAB_MEX_LIBRARY})
    set_target_properties(${${target}_NAME}
        PROPERTIES
        PREFIX ""
        OUTPUT_NAME ${${target}_OUTPUT_NAME}
        SUFFIX "${MEXOCT_MATLAB_EXT}")
    if(WIN32 AND CMAKE_SIZEOF_VOID_P MATCHES "8")
        set(MX_COMPAT "MX_COMPAT_64")
    elseif(WIN32 AND CMAKE_SIZEOF_VOID_P MATCHES "4")
        set (MX_COMPAT "MX_COMPAT_32")
    else()
        set (MX_COMPAT "")
    endif()
    target_compile_definitions(${${target}_NAME} PRIVATE ${MEX_API_MACRO} ${MEX_COMPLEX_INTERLEAVE} ${MX_COMPAT} "USE_MEX_CMD" "MATLAB_MEX_FILE" "MEXOCT_MATLAB=1")

    # entry point in the mex file + taking care of visibility and symbol clashes.
    if(WIN32)
        set_target_properties(${${target}_NAME}
            PROPERTIES
            DEFINE_SYMBOL "DLL_EXPORT_SYM=__declspec(dllexport)")
        if(MINGW)
            set_target_properties(${${target}_NAME} PROPERTIES
                LINK_FLAGS "-static-libgcc -static-libstdc++")
        endif()

    else()

        if(HAS_MINUS_PTHREAD AND NOT APPLE)
            # Apparently, compiling with -pthread generated the proper link flags
            # and some defines at compilation
            target_compile_options(${${target}_NAME} PRIVATE "-pthread")
        endif()


        # if we do not do that, the symbols linked from eg. boost remain weak and
        # then clash with the ones defined in the matlab process. So by default
        # the symbols are hidden.
        # This also means that for shared libraries (like MEX), the entry point
        # should be explicitly declared with default visibility, otherwise Matlab
        # cannot find the entry point.
        # Note that this is particularly meaningful if the MEX wrapper itself
        # contains symbols that are clashing with Matlab (that are compiled in the
        # MEX file). In order to propagate the visibility options to the libraries
        # to which the MEX file is linked against, the -Wl,--exclude-libs,ALL
        # option should also be specified.

        set_target_properties(${${target}_NAME}
            PROPERTIES
            CXX_VISIBILITY_PRESET "hidden"
            C_VISIBILITY_PRESET "hidden"
            VISIBILITY_INLINES_HIDDEN ON
            )


        set_target_properties(${${target}_NAME}
            PROPERTIES
            DEFINE_SYMBOL "DLL_EXPORT_SYM=__attribute__ ((visibility (\"default\")))"
            )

    endif()

endfunction()

