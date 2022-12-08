SET( SRC
    benchmark_all_tgcsylv
    benchmark_all_tglyap
    benchmark_all_tgstein
    benchmark_all_tgsylv
    benchmark_all_trlyap
    benchmark_all_trstein
    benchmark_all_trsylv2
    benchmark_all_trsylv
)

foreach (srcfile ${SRC})
	add_executable (${srcfile} ${srcfile}.c)
    target_link_libraries (
        ${srcfile}
            PRIVATE
                OpenMP::OpenMP_C
                mepack
                ${ADDITIONAL_LIBRARIES}
                benchmark_lib
                ${LAPACK_LIBRARIES}
                ${BLAS_LIBRARIES}
                cscutils)
endforeach (srcfile)



