IF ( CMAKE_Fortran_COMPILER_LOADED )
    INCLUDE(CheckFortranCompilerFlag)

    IF ( CMAKE_Fortran_FLAGS MATCHES "-Wconversion")
        STRING(REPLACE "-Wconversion" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
        SET (Fortran_W_CONVERSION 1)
    ENDIF()

    IF (NOT Fortran_W_CONVERSION )
        check_fortran_compiler_flag("-Wconversion" Fortran_W_CONVERSION)
    ENDIF()

    IF ( Fortran_W_CONVERSION)
        SET (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wconversion")
        SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -Wconversion")
        SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -Wconversion")
        SET (CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} -Wconversion")
    ENDIF()
ENDIF() # CMAKE_Fortran_COMPILER_LOADED


