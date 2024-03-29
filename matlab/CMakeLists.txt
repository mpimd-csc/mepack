SET(MATLAB_MEX_INTERLEAVED_COMPLEX_API ON)

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


FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/matlab")
FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/octave")
FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/help")


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

    IF(WIN32)
        ADD_CUSTOM_COMMAND(TARGET mepackmexoct POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE_NAME:mepackmexoct>" "${PROJECT_BINARY_DIR}/matlab/matlab"
		COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE_NAME:mepackmexoct>" "${PROJECT_BINARY_DIR}/matlab/octave")
    ELSE()
        ADD_CUSTOM_COMMAND(TARGET mepackmexoct POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_SONAME_FILE:mepackmexoct>" "${PROJECT_BINARY_DIR}/matlab/matlab"
            COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_SONAME_FILE:mepackmexoct>" "${PROJECT_BINARY_DIR}/matlab/octave")
    ENDIF()
ENDIF()


SET(MEPACK_MATLAB_TARGETS )
SET(MEPACK_OCTAVE_TARGETS )
SET(MEPACK_OCTAVE_HELP_TARGETS )
SET(MEXOCT_HELP_DIR ${CMAKE_CURRENT_BINARY_DIR}/help)
function (force_static_libs lang)
  set(result )
  foreach (lib ${CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES})
    list(FIND ARGN ${lib} found)
    if (NOT found EQUAL -1)
      find_library(${lib}_LIB
        ${CMAKE_STATIC_LIBRARY_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_${lang}_IMPLICIT_LINK_DIRECTORIES})
      if (${lib}_LIB)
        set(lib ${${lib}_LIB})
      endif ()
    endif ()
    set(result ${result} ${lib})
  endforeach ()
  set(CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES ${result} PARENT_SCOPE)
  message("CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES: ${result}")
endfunction ()


if (WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    force_static_libs(Fortran gfortran quadmath)
    find_library(GOMP NAMES libgomp.a)
else()
    set(GOMP OpenMP::OpenMP_Fortran)
endif()


FOREACH(MO IN LISTS MO_INTERFACES)
    mexoct_add(mepack_${MO}  ${MO}.cpp)
    IF(MEXOCT_LINK_STATIC STREQUAL ON )
        mexoct_target_link_libraries(mepack_${MO} mepack_static ${GOMP})
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
        if (WIN32)
            set_target_properties(mex_mepack_${MO} PROPERTIES LINK_FLAGS "-Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive -Wl,-Bdynamic -static-libgcc -static-libstdc++" )
        endif()
        list(APPEND MEPACK_MATLAB_TARGETS mex_mepack_${MO} )
    endif()

    if (MEXOCT_OCTAVE_FOUND)
        if (NOT MEXOCT_SKIP_RPATH )
            set_target_properties(oct_mepack_${MO} PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                                                      INSTALL_RPATH "$ORIGIN")
        endif()
        list(APPEND MEPACK_OCTAVE_TARGETS oct_mepack_${MO})
        list(APPEND MEPACK_OCTAVE_HELP_TARGETS oct_mepack_${MO}_doc)
     endif()
ENDFOREACH()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_tests_matlab.m.in ${CMAKE_CURRENT_BINARY_DIR}/run_tests_matlab.m @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_tests_octave.m.in ${CMAKE_CURRENT_BINARY_DIR}/run_tests_octave.m @ONLY)

set(MEXOCT_TEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests/")
set(MEXOCT_BUILD_PATH "${CMAKE_CURRENT_BINARY_DIR}")
set(MEXOCT_MATLAB_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/matlab")
set(MEXOCT_OCTAVE_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/octave")

if (WIN32)
    cmake_path(NATIVE_PATH MEXOCT_TEST_PATH NORMALIZE MEXOCT_TEST_PATH)
    cmake_path(NATIVE_PATH MEXOCT_BUILD_PATH NORMALIZE MEXOCT_BUILD_PATH)
    cmake_path(NATIVE_PATH MEXOCT_MATLAB_BINARY_PATH NORMALIZE MEXOCT_MATLAB_BINARY_PATH)
    cmake_path(NATIVE_PATH MEXOCT_OCTAVE_BINARY_PATH NORMALIZE MEXOCT_OCTAVE_BINARY_PATH)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/startup.m.in ${CMAKE_CURRENT_BINARY_DIR}/startup.m @ONLY)

function(COPY_DLL)
    unset(DLL CACHE)
    unset(DLL)
    set(OPTIONS "")
    set(ONE_VALUE "DEST")
    set(MULTI_VALUE "NAMES")
    cmake_parse_arguments(COPY_DLL "${OPTIONS}" "${ONE_VALUE}" "${MULTI_VALUE}" ${ARGN})
    find_file(DLL NAMES ${COPY_DLL_NAMES})
    if (DLL)
        message(STATUS "Add DLL ${DLL} to MATLAB.")
        file(COPY "${DLL}" DESTINATION "${COPY_DLL_DEST}")
    endif()
endfunction()


if (MEXOCT_MATLAB_FOUND)
    if(WIN32)
        add_test(NAME test_matlab COMMAND ${MEXOCT_MATLAB_BINARY} -batch "ierr=run_tests_matlab(); if(ierr); error('Failed'); end;" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    else()
        add_test(NAME test_matlab COMMAND ${MEXOCT_MATLAB_BINARY} -nojvm -nodesktop -nosplash -r "ierr=run_tests_matlab(); quit(ierr);" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    file(GLOB M_HELP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/help/*.m)
    foreach(M_FILE ${M_HELP_FILES})
        file(COPY ${M_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/matlab/)
    endforeach()


    if (WIN32)
        #         if (MINGW)
        #             COPY_DLL(DEST ${CMAKE_CURRENT_BINARY_DIR}/matlab NAMES libgfortran-6.dll libgfortran-5.dll libgfortran-4.dll libgfortran-3.dll)
        #             COPY_DLL(DEST ${CMAKE_CURRENT_BINARY_DIR}/matlab NAMES libgomp-2.dll libgomp-1.dll libgomp.dll)
        #             COPY_DLL(DEST ${CMAKE_CURRENT_BINARY_DIR}/matlab NAMES libwinpthread-2.dll libwinpthread-1.dll libwinpthread.dll)
        #         endif()
        #
        set(WORKDIR_IN "${CMAKE_CURRENT_BINARY_DIR}/matlab" )
        cmake_path(NATIVE_PATH WORKDIR_IN NORMALIZE WORKDIR)
        set(ARCHIVE_IN   "${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-matlab-${MEXOCT_MATLAB_RELEASE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.zip" )
        cmake_path(NATIVE_PATH ARCHIVE_IN NORMALIZE ARCHIVE)

        file(GLOB ARCHIVE_FILES LIST_DIRECTORIES FALSE RELATIVE "${CMAKE_CURRENT_BINARY_DIR}/matlab" "${CMAKE_CURRENT_BINARY_DIR}/matlab/*.m" "${CMAKE_CURRENT_BINARY_DIR}/matlab/*.dll")
        foreach (F ${MO_INTERFACES})
            list(APPEND ARCHIVE_FILES "mepack_${F}${MEXOCT_MATLAB_EXT}")
        endforeach()
        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/dist_files_matlab.lst" "")
        foreach(FN ${ARCHIVE_FILES})
            file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/dist_files_matlab.lst" "${FN}\n")
        endforeach()
        add_custom_target("dist-matlab" ALL
            ${CMAKE_COMMAND} -E tar cf "${ARCHIVE}" --format=zip --files-from="${CMAKE_CURRENT_BINARY_DIR}/dist_files_matlab.lst"
            WORKING_DIRECTORY "${WORKDIR}"
            COMMENT "Create MATLAB distribution ${ARCHIVE}"
            DEPENDS ${MEPACK_MATLAB_TARGETS})
    else()

         # Insert --format=zip to create a zip file
        add_custom_target("dist-matlab" ALL
            ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-matlab-${MEXOCT_MATLAB_RELEASE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz "*"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/matlab"
            COMMENT "Create MATLAB distribution ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-matlab-${MEXOCT_MATLAB_RELEASE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz"
            DEPENDS ${MEPACK_MATLAB_TARGETS})
    endif()
endif()

if (MEXOCT_OCTAVE_FOUND)
    add_test(NAME test_octave COMMAND ${MEXOCT_OCTAVE_EXECUTABLE} --eval "ierr=run_tests_octave(); quit(ierr);" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # Insert --format=zip to create a zip file
    add_custom_target("dist-octave" ALL
        ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-octave-${MEXOCT_OCTAVE_VERSION_STRING}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz "*"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/octave"
        COMMENT "Create OCTAVE distribution ${PROJECT_BINARY_DIR}/mepack-${MEPACK_VERSION}-octave-${MEXOCT_OCTAVE_VERSION_STRING}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-dist.tar.gz"
        DEPENDS ${MEPACK_OCTAVE_TARGETS})

    add_custom_target("update-matlab-doc"
        ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_BINARY_DIR}/help" "${CMAKE_CURRENT_SOURCE_DIR}/help"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Update help texts"
        DEPENDS ${MEPACK_OCTAVE_HELP_TARGETS})

endif()



