set(SRCS
    copy_input
    return_ArgMsg
    square_trait
    type_info
    vector_trait
    vector_trait_column
    vector_trait_row
    )

set(MEXOCT_MATLAB_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/matlab")
set(MEXOCT_OCTAVE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/octave")

configure_file(run_tests_matlab.m.in ${CMAKE_CURRENT_BINARY_DIR}/run_tests_matlab.m @ONLY)
configure_file(run_tests_octave.m.in ${CMAKE_CURRENT_BINARY_DIR}/run_tests_octave.m @ONLY)


foreach(S ${SRCS})
    message(STATUS "Build test ${S}")
    mexoct_add(${S} cpp/${S}.cpp)
endforeach()

if (MEXOCT_MATLAB_FOUND)
    add_test(NAME test_matlab COMMAND ${MEXOCT_MATLAB_BINARY} -nojvm -nodesktop -nosplash -r "ierr=run_tests_matlab(); quit(ierr);" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

if (MEXOCT_OCTAVE_FOUND)
    add_test(NAME test_octave COMMAND ${MEXOCT_OCTAVE_EXECUTABLE} --eval "ierr=run_tests_octave(); quit(ierr);" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

