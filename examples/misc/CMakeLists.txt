SET( SRC
    benchmark_eigenvalue_reorder
    benchmark_dag
    benchmark_row_columnwise
    solve_refine_glyap
    solve_slicot_glyap
    solve_mepack_glyap
)

IF ( RECSY )
    SET( SRC ${SRC} benchmark_recsy)
ENDIF()

foreach (srcfile ${SRC})
	add_executable (${srcfile} ${srcfile}.f90)
    target_link_libraries (
        ${srcfile}
            PRIVATE
            OpenMP::OpenMP_Fortran
                mepack
                ${ADDITIONAL_LIBRARIES}
                benchmark_lib
                ${LAPACK_LIBRARIES}
                ${BLAS_LIBRARIES}
                cscutils)
endforeach (srcfile)



