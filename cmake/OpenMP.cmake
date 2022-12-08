FIND_PACKAGE(OpenMP REQUIRED)

if (CMAKE_C_COMPILER_LOADED)
    if(NOT TARGET OpenMP::OpenMP_C)
        find_package(Threads REQUIRED)
        add_library(OpenMP::OpenMP_C IMPORTED INTERFACE)
        set_property(TARGET OpenMP::OpenMP_C
            PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_C_FLAGS})
        # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
        set_property(TARGET OpenMP::OpenMP_C
            PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_C_FLAGS} Threads::Threads)

    endif()
endif()

if ( CMAKE_CXX_COMPILER_LOADED )
    # For CMake < 3.9, we need to make the target ourselves
    if(NOT TARGET OpenMP::OpenMP_CXX)
        find_package(Threads REQUIRED)
        add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
        set_property(TARGET OpenMP::OpenMP_CXX
            PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
        # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
        set_property(TARGET OpenMP::OpenMP_CXX
            PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

    endif()
endif()

if ( CMAKE_Fortran_COMPILER_LOADED)
    if(NOT TARGET OpenMP::OpenMP_Fortran)
        find_package(Threads REQUIRED)
        add_library(OpenMP::OpenMP_Fortran IMPORTED INTERFACE)
        set_property(TARGET OpenMP::OpenMP_Fortran
            PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_Fortran_FLAGS})
        # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
        set_property(TARGET OpenMP::OpenMP_Fortran
            PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_Fortran_FLAGS} Threads::Threads)

    endif()
endif()

