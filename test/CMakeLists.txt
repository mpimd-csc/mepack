ADD_SUBDIRECTORY(arguments)


add_executable(memory memory.c)
target_link_libraries (memory mepack)

add_test(MEPACK_MEMORY ${CMAKE_CURRENT_BINARY_DIR}/memory)

