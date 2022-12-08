if(NOT DEFINED CMAKE_INSTALL_MODULEDIR)
    if(NOT DEFINED Fortran_MODULE_NAME)
        set(Fortran_MODULE_NAME ${PROJECT_NAME})
    endif()
  set(
    CMAKE_INSTALL_MODULEDIR
    "${CMAKE_INSTALL_INCLUDEDIR}/${Fortran_MODULE_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
    CACHE
    STRING
    "Directory in prefix to install generated module files"
  )
endif()

