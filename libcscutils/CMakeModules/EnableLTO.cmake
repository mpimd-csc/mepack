IF ( CMAKE_C_COMPILER_LOADED )
    INCLUDE(CheckCCompilerFlag)

    IF ( CMAKE_C_FLAGS MATCHES "-flto=auto")
        STRING(REPLACE "-flto=auto" "" CMAKE_C_FLAGS_X "${CMAKE_C_FLAGS}")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_X}" CACHE INTERNAL "")
        SET (C_F_LTO 1)
    ENDIF()

    IF (NOT C_F_LTO )
        check_c_compiler_flag("-flto=auto" C_F_LTO)
    ENDIF()

    IF ( C_F_LTO)
        SET (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -flto=auto" CACHE INTERNAL "")
    ENDIF()
ENDIF() # CMAKE_C_COMPILER_LOADED

IF ( CMAKE_CXX_COMPILER_LOADED )
    INCLUDE(CheckCXXCompilerFlag)

    IF ( CMAKE_CXX_FLAGS MATCHES "-flto=auto")
        STRING(REPLACE "-flto=auto" "" CMAKE_CXX_FLAGS_X "${CMAKE_CXX_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_X}" CACHE INTERNAL "")

        SET (CXX_F_LTO 1)
    ENDIF()

    IF (NOT CXX_F_LTO )
        check_cxx_compiler_flag("-flto=auto" CXX_F_LTO)
    ENDIF()

    IF ( CXX_F_LTO)
        SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -flto=auto" CACHE INTERNAL "")
    ENDIF()
ENDIF() # CMAKE_CXX_COMPILER_LOADED

IF ( CMAKE_Fortran_COMPILER_LOADED )
    INCLUDE(CheckFortranCompilerFlag)

    IF ( CMAKE_Fortran_FLAGS MATCHES "-flto=auto")
        STRING(REPLACE "-flto=auto" "" CMAKE_Fortran_FLAGS_X "${CMAKE_Fortran_FLAGS}")
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_X}" CACHE INTERNAL "")

        SET (Fortran_F_LTO 1)
    ENDIF()

    IF (NOT Fortran_F_LTO )
        check_fortran_compiler_flag("-flto=auto" Fortran_F_LTO)
    ENDIF()

    IF ( Fortran_F_LTO)
        SET (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -flto=auto" CACHE INTERNAL "")
        SET (CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} -flto=auto" CACHE INTERNAL "")
    ENDIF()
ENDIF() # CMAKE_Fortran_COMPILER_LOADED


