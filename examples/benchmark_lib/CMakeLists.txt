SET(SRC init.c
	problems.c
	rhs.c
	recsy_helper.c
)


FIND_PACKAGE(HDF5 REQUIRED)

ADD_SUBDIRECTORY(slicot)

ADD_LIBRARY(benchmark_lib STATIC ${SRC} $<TARGET_OBJECTS:slicot> )
TARGET_INCLUDE_DIRECTORIES(
    benchmark_lib
        PRIVATE
            ${HDF5_INCLUDE_DIRS}
        PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/examples/benchmark_lib>
    )
TARGET_LINK_LIBRARIES(
    benchmark_lib
        PUBLIC
            cscutils
            mepack
    )

