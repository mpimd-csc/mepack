
set(MEXOCT_MATLAB_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/matlab")
set(MEXOCT_OCTAVE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/octave")

mexoct_add(example_info example_info.cpp)
mexoct_add(allocate_memory allocate_memory.cpp)
mexoct_add(no_argument no_argument.cpp)
mexoct_add(example_function_variants example_function_variants.cpp)



#
# Build the QR example
#
mexoct_add(geqrf geqrf.cpp)
mexoct_target_link_blas(geqrf)

