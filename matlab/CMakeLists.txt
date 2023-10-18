INCLUDE(${PROJECT_SOURCE_DIR}/mexoct/MexOct.cmake)

IF(MEXOCT_MATLAB_FOUND STREQUAL ON AND INTEGER8 STREQUAL OFF)
    MESSAGE(FATAL_ERROR "If the MATLAB interface should be built, 64 bit integers have to used. Enable them by specifying -DINTEGER8=ON or disable building the MATLAB interface by -DMEXOCT_MATLAB=OFF")
ENDIF()
IF(MEXOCT_OCTAVE_FOUND AND INTEGER8 STREQUAL ON )
    MESSAGE(WARNING "MEPACK is built with 64 intergers but Octave's BLAS and LAPACK libraries might not. Please check this or disable the 64 bit intergers using -DINTEGER8=OFF or disable the Octave interface with -DMEXOCT_OCTAVE=OFF")
ENDIF()

SET(MO_INTERFACES
    trlyap
    trstein
    tglyap
    tgstein
    trsylv
    trsylv2
    tgsylv
    tgcsylv
    tgcsylv_dual
    lyap
    stein
    glyap
    gstein
    sylv
    sylv2
    gsylv
    csylv
    csylv_dual
    lyap_refine
    stein_refine
    glyap_refine
    gstein_refine
    sylv_refine
    sylv2_refine
    gsylv_refine
    csylv_refine
    csylv_dual_refine
    )

OPTION(MEXOCT_LINK_STATIC "Static Linking of MEX/OCT files" OFF)
IF (MEXOCT_LINK_STATIC STREQUAL ON AND BUILD_STATIC STREQUAL OFF)
    MESSAGE(FATAL_ERROR "The MATLAB/Octave Interface can only be build with static libraries of BUILD_STATIC=ON")
ENDIF()

IF ( NOT DEFINED MEXOCT_SKIP_RPATH )
    SET(MEXOCT_SKIP_RPATH false)
ENDIF()

# Link Library
IF (MEXOCT_LINK_STATIC STREQUAL OFF)
    ADD_LIBRARY(mepackmexoct SHARED $<TARGET_OBJECTS:mepack> $<TARGET_OBJECTS:lapack_missing> $<TARGET_OBJECTS:cscutils_lua_src> $<TARGET_OBJECTS:lua_src>)
    TARGET_LINK_LIBRARIES(
        mepackmexoct
        PUBLIC
        OpenMP::OpenMP_Fortran
        )
    TARGET_INCLUDE_DIRECTORIES(
        mepackmexoct
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src/>
        )

    #
    # Default LUA Config
    #
    # INCLUDE(${PROJECT_SOURCE_DIR}/libcscutils/src/lua/LuaInterface.cmake)
    # IF ( NOT MEPACK_LUA_CONFIG )
        # CSC_LUA_ADD_INLINE_CODE(mepackmexoct ${PROJECT_SOURCE_DIR}/src/default_config.lua VARNAME mepack_default_config_code FUNCTIONNAME mepack_default_config)
    # ELSE()
        # CSC_LUA_ADD_INLINE_CODE(mepackmexoct ${MEPACK_LUA_CONFIG} VARNAME mepack_default_config_code FUNCTIONNAME mepack_default_config)
    # ENDIF()

    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/matlab")
    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/octave")

    ADD_CUSTOM_COMMAND(TARGET mepackmexoct POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_SONAME_FILE:mepackmexoct>" "${PROJECT_BINARY_DIR}/matlab/matlab"
        COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_SONAME_FILE:mepackmexoct>" "${PROJECT_BINARY_DIR}/matlab/octave")

ENDIF()


SET(MEPACK_MATLAB_TARGETS )
SET(MEPACK_OCTAVE_TARGETS )

FOREACH(MO IN LISTS MO_INTERFACES)
    mexoct_add(mepack_${MO}  ${MO}.cpp)
    IF(MEXOCT_LINK_STATIC STREQUAL ON )
        mexoct_target_link_libraries(mepack_${MO} mepack_static)
    ELSE()
        mexoct_target_link_libraries(mepack_${MO} mepackmexoct)
    ENDIF()

    mexoct_target_link_blas(mepack_${MO})
    mexoct_target_include_directories(mepack_${MO} PRIVATE ${PROJECT_SOURCE_DIR}/src/)

    if (MEXOCT_MATLAB_FOUND)
        if (NOT MEXOCT_SKIP_RPATH )
            set_target_properties(mex_mepack_${MO} PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                INSTALL_RPATH "$ORIGIN")
        endif()
        list(APPEND MEPACK_MATLAB_TARGETS mex_mepack_${MO} mex_mepack_${MO}_doc)
    endif()

    if (MEXOCT_OCTAVE_FOUND)
        if (NOT MEXOCT_SKIP_RPATH )
            set_target_properties(oct_mepack_${MO} PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                                                      INSTALL_RPATH "$ORIGIN")
        endif()
        list(APPEND MEPACK_OCTAVE_TARGETS oct_mepack_${MO})
     endif()
ENDFOREACH()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_tests_matlab.m.in ${CMAKE_CURRENT_BINARY_DIR}/run_tests_matlab.m @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_tests_octave.m.in ${CMAKE_CURRENT_BINARY_DIR}/run_tests_octave.m @ONLY)

set(MEXOCT_TEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests/")
set(MEXOCT_BUILD_PATH "${CMAKE_CURRENT_BINARY_DIR}")
set(MEXOCT_MATLAB_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/matlab")
set(MEXOCT_OCTAVE_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/octave")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/startup.m.in ${CMAKE_CURRENT_BINARY_DIR}/startup.m @ONLY)

if (MEXOCT_MATLAB_FOUND)
    add_test(NAME test_matlab COMMAND ${MEXOCT_MATLAB_BINARY} -nojvm -nodesktop -nosplash -r "ierr=run_tests_matlab(); quit(ierr);" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # Insert --format=zip to create a zip file
    add_custom_target("dist-matlab" ALL
        ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-matlab-${MEXOCT_MATLAB_RELEASE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz "*"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/matlab"
        COMMENT "Create MATLAB distribution ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-matlab-${MEXOCT_MATLAB_RELEASE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz"
        DEPENDS ${MEPACK_MATLAB_TARGETS})

endif()

if (MEXOCT_OCTAVE_FOUND)
    add_test(NAME test_octave COMMAND ${MEXOCT_OCTAVE_EXECUTABLE} --eval "ierr=run_tests_octave(); quit(ierr);" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # Insert --format=zip to create a zip file
    add_custom_target("dist-octave" ALL
        ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-octave-${MEXOCT_OCTAVE_VERSION_STRING}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz "*"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/octave"
        COMMENT "Create OCTAVE distribution ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-octave-${MEXOCT_OCTAVE_VERSION_STRING}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz"
        DEPENDS ${MEPACK_OCTAVE_TARGETS})
endif()



