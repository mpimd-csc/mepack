SET (SRC
	generate_gevp
	generate_evp
)

foreach (srcfile ${SRC})
	add_executable (${srcfile} ${srcfile}.c)
    target_link_libraries (${srcfile}
        PRIVATE
            OpenMP::OpenMP_C
            OpenMP::OpenMP_Fortran
            mepack
            ${ADDITIONAL_LIBRARIES}
            benchmark_lib
            cscutils
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES}
    )
endforeach (srcfile)
