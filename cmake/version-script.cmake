include(CheckCCompilerFlag)
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18.0)
  include(CheckLinkerFlag)
endif()

macro(_CHECK_VERSION_SCRIPT_SUPPORT)
    set(VERSION_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/version-script-test.lds")
    set(VERSION_SCRIPT_CONTENT "TEST_1.0 {\n}\;")
    file(WRITE "${VERSION_SCRIPT}" ${VERSION_SCRIPT_CONTENT})
    if(COMMAND check_linker_flag)
        check_linker_flag("C" "-Wl,--version-script=${VERSION_SCRIPT}"
            HAS_VERSION_SCRIPT_SUPPORT)
    else()
        check_c_compiler_flag("-Wl,--version-script=${VERSION_SCRIPT}"
            HAS_VERSION_SCRIPT_SUPPORT)
    endif()
    set(_HAS_VERSION_SCRIPT_SUPPORT
        ${HAS_VERSION_SCRIPT_SUPPORT}
        CACHE INTERNAL "Linker supports version scripts")
    # file(REMOVE "${VERSION_SCRIPT}")
endmacro(_CHECK_VERSION_SCRIPT_SUPPORT)

macro(ADD_VERSION_SCRIPT TARGET VERSION_SCRIPT)
  if(NOT DEFINED _HAS_VERSION_SCRIPT_SUPPORT)
    _check_version_script_support()
  endif()
  if(_HAS_VERSION_SCRIPT_SUPPORT)
    if(TARGET ${TARGET})
      set_property(
        TARGET ${TARGET}
        APPEND_STRING
        PROPERTY LINK_FLAGS " -Wl,--version-script=${VERSION_SCRIPT}")
      set_target_properties(${TARGET} PROPERTIES LINK_DEPENDS ${VERSION_SCRIPT})
    endif()
  endif()
endmacro(ADD_VERSION_SCRIPT)

