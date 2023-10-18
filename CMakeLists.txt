CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libcscutils/CMakeModules")
INCLUDE(CheckPolicy)



# Set CMAKE Policies
CheckPolicy(CMP0048 NEW)
CheckPolicy(CMP0069 NEW)
CheckPolicy(CMP0075 NEW)
CheckPolicy(CMP0115 NEW)

# General Project Setup
SET(VERSION 1.1.0)
PROJECT(MEPACK
        VERSION ${VERSION}
	    LANGUAGES C Fortran
        DESCRIPTION "MEPACK is an Fortran software library for the solution of dense Sylvester-likematrix equations.")

IF ( NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    MESSAGE(STATUS "MEPACK: Build within another project. Tests, examples, documentation, and the MATLAB interface are disabled.")
    SET(MEPACK_IN_PROJECT 1)
ELSE()
    SET(MEPACK_IN_PROJECT 0)
ENDIF()

INCLUDE(GNUInstallDirs)

# Options
OPTION(DEBUG "Enable Debug" OFF)
OPTION(LTO  "Enable Link Time Optimization" ON)
OPTION(HOSTOPT "Enable Host Optimization" OFF)
# OPTION(BUILD_SHARED_LIBS "Build a shared libary" ON)
OPTION(INTEGER8 "Build 8-byte integer aware libary." OFF)
OPTION(EXAMPLES "Build Examples" ON)
OPTION(MATLAB "Build the MATLAB and Octave Interface" OFF)
OPTION(DOC "Build documentation" OFF)
OPTION(BUILD_STATIC "Build Static Libraries" OFF)

SET(CMAKE_POSITION_INDEPENDENT_CODE True)


IF ( DEBUG STREQUAL ON )
	SET (CMAKE_BUILD_TYPE "Debug")
ELSEIF(NOT CMAKE_BUILD_TYPE)
	SET (CMAKE_BUILD_TYPE "Release")
ENDIF()

IF (INTEGER8 STREQUAL ON)
    SET(MEPACK_SUFFIX 64)
    SET(MEPACK_NAME mepack64)
ELSE()
    SET(MEPACK_SUFFIX "")
    SET(MEPACK_NAME mepack)
ENDIF()

IF (MATLAB STREQUAL ON)
    ENABLE_LANGUAGE(CXX)
ENDIF()

# Find Valgrind for testing
IF ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    FIND_PROGRAM(VALGRIND_TOOL valgrind)
    IF ( VALGRIND_TOOL )
        MESSAGE(STATUS "Found valgrind: ${VALGRIND_TOOL}")
        SET(MEMORYCHECK_COMMAND "${VALGRIND_TOOL}")
        SET(CTEST_MEMORYCHECK_COMMAND "${VALGRIND_TOOL}")
    ENDIF()
    IF ( VALGRIND_SUPPRESSIONS )
        IF(NOT IS_ABSOLUTE ${VALGRIND_SUPPRESSIONS})
            SET(VALGRIND_SUPPRESSIONS "${CMAKE_CURRENT_BINARY_DIR}/${VALGRIND_SUPPRESSIONS}")
        ENDIF()
        MESSAGE(STATUS "Use valgrind suppressions from: ${VALGRIND_SUPPRESSIONS}")
        SET(MEMORYCHECK_SUPPRESSIONS_FILE "${VALGRIND_SUPPRESSIONS}" CACHE INTERNAL "X")
        SET(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "${VALGRIND_SUPPRESSIONS}" CACHE INTERNAL "X")
    ELSE()
        MESSAGE(STATUS "No suppression for valgrind given. This leads to false positives during the memory check.")
    ENDIF()
ENDIF()

ENABLE_TESTING()
INCLUDE(CTest)

# Add OpenMP support
INCLUDE(cmake/OpenMP.cmake)

# Fortran Module Support
SET(Fortran_MODULE_NAME ${MEPACK_NAME})
INCLUDE(cmake/FortranModuleDir.cmake)
INCLUDE(cmake/version-script.cmake)

IF (OPENMP_FOUND)
	MESSAGE(STATUS "OpenMP : ${OpenMP_Fortran_VERSION}")
    SET(OpenMP_Fortran_VERSION ${OpenMP_Fortran_VERSION} CACHE INTERNAL "OpenMP Fortran Version")
    IF ("${OpenMP_Fortran_VERSION}" VERSION_LESS "4.0")
        MESSAGE(WARNING "The Fortran Compiler does not fully upport OpenMP 4 or the version could not be detected correctly. This may cause problems.")
    ENDIF()
ENDIF()

INCLUDE(FortranCompilerSettings)
INCLUDE(CCompilerSettings)
# INCLUDE(CXXCompilerSettings)
INCLUDE(CMakeHelpers)

# Name Mangling
include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
	                 MACRO_NAMESPACE "FC_"
			 SYMBOL_NAMESPACE "FC_")
FortranCInterface_VERIFY()

# BLAS and LAPACK
IF ( INTEGER8 STREQUAL ON)
	FIND_PACKAGE(BLAS64 REQUIRED)
    FIND_PACKAGE(LAPACK64 REQUIRED)
    SET (BLAS_LIBRARIES ${BLAS64_LIBRARIES})
    SET (LAPACK_LIBRARIES ${LAPACK64_LIBRARIES})
    MESSAGE(STATUS "Build with ILP64 memory model")
    SET(BLAS64_FOUND true)
    SET(LAPACK64_FOUND true)
ELSE()
	FIND_PACKAGE(BLAS REQUIRED)
	FIND_PACKAGE(LAPACK REQUIRED)
ENDIF()
SET(LIBRARIES ${LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})



# LAPACK setup
SET(LAPACK_MISSING_TARGET_NAME lapack_missing)
INCLUDE(${PROJECT_SOURCE_DIR}/lapack-missing/LapackMissing.cmake)
# SET(LIBRARIES lapack_missing ${LIBRARIES})

CHECK_ADD_C_FLAG("-Wno-overlength-strings" NO_OVERLENGTH)

# Make the compilers working more strict
INCLUDE(EnablePedanticC)
INCLUDE(EnablePedanticFortran)
INCLUDE(EnableStrictPrototypesC)
INCLUDE(EnableCastQual)
INCLUDE(EnableMayBeUninitialized)
INCLUDE(EnableStringOpTruncation)
INCLUDE(EnableConversionFortran)
IF (LTO STREQUAL ON)
    INCLUDE(EnableLTO)
ENDIF()


# COVERAGE TEST
IF(COVERAGE)
    # check if we use a gnu compiler
    IF(NOT CMAKE_C_COMPILER_ID  STREQUAL "GNU")
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Code coverage is only possible with GNU Compilers")
    ENDIF()

    # add flags
    SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} gcov)
    CHECK_ADD_C_FLAG("-g" SUCCESS1)
    SET(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -g")
    CHECK_ADD_C_FLAG("-ftest-coverage" SUCCESS3)
    SET(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -ftest-coverage")
    CHECK_ADD_C_FLAG("-fprofile-arcs" SUCCESS2)
    SET(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs")

    CHECK_ADD_Fortran_FLAG("-g" SUCCESS4)
    CHECK_ADD_Fortran_FLAG("-ftest-coverage" SUCCESS6)
    CHECK_ADD_Fortran_FLAG("-fprofile-arcs" SUCCESS5)
    SET(CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -g -fprofile-arcs -ftest-coverage")

    IF(NOT(SUCCESS1 AND SUCCESS2 AND SUCCESS3 AND SUCCESS4 AND SUCCESS5 AND SUCCESS6))
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Compiler Flags for coverage report are not supported by your compiler")
    ENDIF()

    # find gcov library, gcov lcov and genhtml program and add new target make coverage
    FIND_PROGRAM(LCOV       lcov)
    FIND_PROGRAM(GENHTML    genhtml)
    FIND_PROGRAM(GCOV       gcov)
    FIND_PROGRAM(SED        sed)
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found lcov: ${LCOV}")
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found gcov: ${GCOV}")
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found genhtml: ${GENHTML}")
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found sed: ${SED}")


    IF(LCOV AND GENHTML AND GCOV AND SED)
        # add target make coverage
        ADD_CUSTOM_TARGET(
            coverage
            COMMAND ${LCOV} --zerocounters --directory "${CMAKE_BINARY_DIR}/"
            COMMAND make -C "${CMAKE_BINARY_DIR}" test
            COMMAND ${LCOV} -c --gcov-tool ${GCOV} --output-file "${CMAKE_BINARY_DIR}/lcov.info" --directory "${CMAKE_BINARY_DIR}/"
            COMMAND ${SED} -i -e 's,${CMAKE_BINARY_DIR}/src/CMakeFiles/mepack.dir/,${PROJECT_SOURCE_DIR}/src/,g' -e 's,${CMAKE_BINARY_DIR}/libcscutils/src/quad.f90,${PROJECT_SOURCE_DIR}/libcscutils/src/common/quad.f90,g' "${CMAKE_BINARY_DIR}/lcov.info"
            COMMAND ${GENHTML} --legend --output "${CMAKE_BINARY_DIR}/coverage" "${CMAKE_BINARY_DIR}/lcov.info"
            COMMAND  ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Coverage Report available: ${CMAKE_BINARY_DIR}/coverage"
            )
        ADD_CUSTOM_TARGET(
            coverage-norun
            COMMAND ${LCOV} -c --gcov-tool ${GCOV} --output-file "${CMAKE_BINARY_DIR}/lcov.info" --directory "${CMAKE_BINARY_DIR}/"
            COMMAND ${SED} -i -e 's,${CMAKE_BINARY_DIR}/src/CMakeFiles/mepack.dir/,${PROJECT_SOURCE_DIR}/src/,g' -e 's,${CMAKE_BINARY_DIR}/libcscutils/src/quad.f90,${PROJECT_SOURCE_DIR}/libcscutils/src/common/quad.f90,g' "${CMAKE_BINARY_DIR}/lcov.info"
            COMMAND ${GENHTML} --legend --output "${CMAKE_BINARY_DIR}/coverage" "${CMAKE_BINARY_DIR}/lcov.info"
            COMMAND  ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Coverage Report available: ${CMAKE_BINARY_DIR}/coverage"
            )

        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make coverage")
        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make coverage-norun")
    ELSE()
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Coverage will not work due to missing tools. Need lcov, gcov and genhtml.")
    ENDIF()

ENDIF()

IF (EXAMPLES AND NOT MEPACK_IN_PROJECT)
    SET(CSCUTILS_FEATURES "hdf5;hardware;io;inifile;lua")
ELSE()
    SET(CSCUTILS_FEATURES "lua")
ENDIF()



ADD_SUBDIRECTORY(libcscutils)

ADD_SUBDIRECTORY(src)
IF ( EXAMPLES AND NOT MEPACK_IN_PROJECT)
    ADD_SUBDIRECTORY(examples)
    ADD_SUBDIRECTORY(test)
ENDIF()
IF(MATLAB AND NOT MEPACK_IN_PROJECT)
    ADD_SUBDIRECTORY(matlab)
ENDIF()

# Documentation
IF (DOC AND NOT MEPACK_IN_PROJECT)
    FIND_PACKAGE(Doxygen REQUIRED dot)
    FIND_PROGRAM(SED_TOOL NAMES sed gsed)
    IF ( NOT SED_TOOL )
        MESSAGE(FATAL_ERROR "Building the documentation requires sed / gsed")
    ENDIF()

    IF (DOXYGEN_FOUND)
        IF (DOXYGEN_VERSION VERSION_LESS "1.8.16")
            MESSAGE(FATAL_ERROR "In order to create the documentation doxygen >= 1.8.16 is required. Please update doxygen or turn of the documentation with -DDOC=OFF")
        ENDIF()

        IF (DOXYGEN_VERSION VERSION_EQUAL "1.8.17")
            MESSAGE(STATUS "Doxygen ${DOXYGEN_VERSION} does not play well with Ghostscript. MATHJAX is used to display formulas in the documentation.")
            SET(USE_MATHJAX "YES")
            MESSAGE(STATUS "Download MathJax 2.7.x is required due to compatability issues.")
            FIND_PROGRAM(WGET NAMES wget)
            FIND_PROGRAM(TAR NAMES tar gtar)
            IF ( NOT WGET )
                MESSAGE(FATAL_ERROR "wget is required to download MathJax")
            ENDIF ()
            IF ( NOT TAR )
                MESSAGE(FATAL_ERROR "tar is required to extract MathJax")
            ENDIF ()
            # FILE(DOWNLOAD https://github.com/mathjax/MathJax/archive/refs/tags/2.7.9.tar.gz ${CMAKE_BINARY_DIR}/doc/mathjax.tar.gz SHOW_PROGRESS)
            FILE(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
            ADD_CUSTOM_COMMAND(OUTPUT "${CMAKE_BINARY_DIR}/doc/explore-html/mathjax/MathJax.js"
                COMMAND ${WGET} -c https://github.com/mathjax/MathJax/archive/refs/tags/2.7.9.tar.gz -O ${CMAKE_BINARY_DIR}/doc/mathjax.tar.gz
                COMMAND ${TAR} -x --transform="s,MathJax-2.7.9,explore-html/mathjax,g" -f "${CMAKE_BINARY_DIR}/doc/mathjax.tar.gz"
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc
                )
            ADD_CUSTOM_TARGET(create_mathjax DEPENDS "${CMAKE_BINARY_DIR}/doc/explore-html/mathjax/MathJax.js")
        ELSE()
            SET(USE_MATHJAX "NO")
            ADD_CUSTOM_TARGET(create_mathjax)
        ENDIF()

        # set input and output files
        set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        FILE(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/matlab")
        add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/README.post.md"
            COMMAND ${SED_TOOL} -e 's/\\$$`/\\\\f$$/g' -e 's/`\\$$/\\\\f$$/g' -e 's,doc/install.md,@ref install,g' -e 's,doc/tuning.md,@ref tuning,g'  "${PROJECT_SOURCE_DIR}/README.md" > "${CMAKE_BINARY_DIR}/README.post.md"
            DEPENDS "${PROJECT_SOURCE_DIR}/README.md"
            )
        add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/matlab/README.post.md"
            COMMAND ${SED_TOOL} -e 's/\\$$`/\\\\f$$/g' -e 's/`\\$$/\\\\f$$/g' -e 's,../doc/install.md,@ref install,g' -e 's,doc/tuning.md,@ref tuning,g'  "${PROJECT_SOURCE_DIR}/matlab/README.md" > "${CMAKE_BINARY_DIR}/matlab/README.post.md"
            DEPENDS "${PROJECT_SOURCE_DIR}/matlab/README.md"
            )

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            DEPENDS create_mathjax "${CMAKE_BINARY_DIR}/README.post.md" "${CMAKE_BINARY_DIR}/matlab/README.post.md"
            VERBATIM )
    ENDIF()
ENDIF()

# Install
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    MEPACK${MEPACK_SUFFIX}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)
install(EXPORT MEPACKTargets
    FILE MEPACK${MEPACK_SUFFIX}Targets.cmake
    NAMESPACE MEPACK::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mepack${MEPACK_SUFFIX}
)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/MEPACKConfig.cmake.in MEPACK${MEPACK_SUFFIX}Config.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mepack${MEPACK_SUFFIX})
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/MEPACK${MEPACK_SUFFIX}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/MEPACK${MEPACK_SUFFIX}ConfigVersion.cmake"
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mepack${MEPACK_SUFFIX}
        )
# PKGconfig

# Export a pkg-config file
IF (INTEGER8 STREQUAL ON )
    SET(PKG_CONFIG_ADD_FLAGS "${I8FLAG}")
ELSE()
    SET(PKG_CONFIG_ADD_FLAGS "")
ENDIF()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mepack.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${MEPACK_NAME}.pc"
  @ONLY
)
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${MEPACK_NAME}.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)


# Final

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" CACHE STRING "Fortran Flags")
SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" CACHE STRING "Fortran Flags DEBUG")
SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "Fortran Flags RELEASE")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "C Flags")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING "C Flags DEBUG")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "C Flags RELEASE")


MESSAGE(STATUS "--------------------")
MESSAGE(STATUS "Final Configuration:")
MESSAGE(STATUS "CMAKE_BUILD_TYPE:            ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS:         ${CMAKE_Fortran_FLAGS}")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS_DEBUG:   ${CMAKE_Fortran_FLAGS_DEBUG}")
MESSAGE(STATUS "CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
MESSAGE(STATUS "CMAKE_C_FLAGS:               ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "CMAKE_C_FLAGS_DEBUG:         ${CMAKE_C_FLAGS_DEBUG}")
MESSAGE(STATUS "CMAKE_C_FLAGS_RELEASE:       ${CMAKE_C_FLAGS_RELEASE}")
MESSAGE(STATUS "BLAS_LIBRARIES:              ${BLAS_LIBRARIES}")
MESSAGE(STATUS "LAPACK_LIBRARIES:            ${LAPACK_LIBRARIES}")
MESSAGE(STATUS "INTEGER8:                    ${INTEGER8}")
MESSAGE(STATUS "MATLAB/Octave interface:     ${MATLAB}")
MESSAGE(STATUS "--------------------")


