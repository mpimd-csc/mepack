IF ( CMAKE_CXX_COMPILER_LOADED )
    INCLUDE(CheckCXXCompilerFlag)
    IF ( CMAKE_CXX_FLAGS MATCHES "-Wpedantic")
        STRING(REPLACE "-Wpedantic" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        SET (CXX_W_PEDANTIC 1)
    ENDIF()

    IF (NOT CXX_W_PEDANTIC)
        check_cxx_compiler_flag("-Wpedantic" CXX_W_PEDANTIC)
    ENDIF()

    IF (CXX_W_PEDANTIC)
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wpedantic")
        SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wpedantic")
        SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wpedantic")
        SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Wpedantic")
    ENDIF()
ENDIF() # CMAKE_CXX_COMPILER_LOADED


