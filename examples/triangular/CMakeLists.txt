SET( SRC

    benchmark_trsylv
	benchmark_trsylv2
	benchmark_trlyap
	benchmark_trstein
	benchmark_tgsylv
	benchmark_tglyap
	benchmark_tgstein
	benchmark_tgcsylv
	benchmark_tgcsylv_dual


	benchmark_trlyap_single
	benchmark_trsylv_single
	benchmark_trsylv2_single
	benchmark_trstein_single
	benchmark_tgsylv_single
	benchmark_tglyap_single
	benchmark_tgstein_single
	benchmark_tgcsylv_single
	benchmark_tgcsylv_dual_single


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


# Testing
# SET(ROWS 1027)
# SET(COLS 1457)
SET(ROWS 257)
SET(COLS 263)
SET(MB 64)
SET(NB 32)
SET(NMAT 1)
SET(ROWS_BIG 527)
SET(ROWS_SMALL 27)
SET(COLS_BIG 607)
SET(COLS_SMALL 21)
#
# Test with M=N matrices for LYAP and STEIN
#
# TRLYAP
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24 )
	ADD_TEST(DOUBLE_TRLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TRLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24)
	ADD_TEST(SINGLE_TRLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

# Big
FOREACH(SOLVER 13 14 15 16 17 18)
	ADD_TEST(DOUBLE_TRLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TRLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trlyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

#
# TGLYAP
#
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24  27 )
	ADD_TEST(DOUBLE_TGLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)

ENDFOREACH()
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24 )
	ADD_TEST(SINGLE_TGLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 13 14 15 16 17 18)
	ADD_TEST(DOUBLE_TGLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGLYAP_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tglyap_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()




# TRSTEIN and TGSTEIN
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24  )
	ADD_TEST(DOUBLE_TRSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TRSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24)
	ADD_TEST(SINGLE_TRSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single    --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TRSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single    --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()

FOREACH(SOLVER 13 14 15 16 17 18)
	ADD_TEST(DOUBLE_TRSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TRSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TRSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)

	ADD_TEST(SINGLE_TRSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single    --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS_BIG}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TRSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single    --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TRSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trstein_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGSTEIN_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgstein_single   --nmat=${NMAT} --rows=${ROWS_SMALL}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()


#
# Test with M!=N for Tests with Sylv
#
# TRSYLV (SIGN = 1)
FOREACH(SOLVER RANGE 1 24 1)
    IF (SOLVER LESS 13 OR SOLVER GREATER 18 )
        SET(ROWSX ${ROWS})
        SET(COLSX ${COLS})
    ELSE()
        SET(ROWSX ${ROWS_BIG})
        SET(COLSX ${COLS_BIG})
    ENDIF()

    ADD_TEST(DOUBLE_TRSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)

    ADD_TEST(DOUBLE_TRSYLV_N_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
	ADD_TEST(DOUBLE_TRSYLV_N_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1 )
	ADD_TEST(DOUBLE_TRSYLV_T_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
	ADD_TEST(DOUBLE_TRSYLV_T_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1 )


    # SMALL
    ADD_TEST(DOUBLE_TRSYLV_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)

    ADD_TEST(DOUBLE_TRSYLV_N_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
    ADD_TEST(DOUBLE_TRSYLV_N_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1 )
    ADD_TEST(DOUBLE_TRSYLV_T_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
    ADD_TEST(DOUBLE_TRSYLV_T_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1 )


    # SYLV2
    ADD_TEST(DOUBLE_TRSYLV2_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV2_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV2_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV2_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TRSYLV2_N_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TRSYLV2_N_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TRSYLV2_T_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TRSYLV2_T_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)

    ADD_TEST(DOUBLE_TRSYLV2_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV2_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV2_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV2_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(DOUBLE_TRSYLV2_N_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(DOUBLE_TRSYLV2_N_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(DOUBLE_TRSYLV2_T_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(DOUBLE_TRSYLV2_T_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2 --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)


ENDFOREACH()


FOREACH(SOLVER RANGE 1 25 1)
	IF (SOLVER LESS 13 OR SOLVER GREATER 18 )
        SET(ROWSX ${ROWS})
        SET(COLSX ${COLS})
    ELSE()
        SET(ROWSX ${ROWS_BIG})
        SET(COLSX ${COLS_BIG})
    ENDIF()
    ADD_TEST(SINGLE_TRSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV_N_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
	ADD_TEST(SINGLE_TRSYLV_N_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1 )
	ADD_TEST(SINGLE_TRSYLV_T_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
	ADD_TEST(SINGLE_TRSYLV_T_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1 )

    #SMALL
    ADD_TEST(SINGLE_TRSYLV_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV_N_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
    ADD_TEST(SINGLE_TRSYLV_N_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1 )
    ADD_TEST(SINGLE_TRSYLV_T_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1 )
    ADD_TEST(SINGLE_TRSYLV_T_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1 )


    IF (NOT SOLVER EQUAL 25)
    ADD_TEST(SINGLE_TRSYLV2_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV2_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV2_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV2_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TRSYLV2_N_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TRSYLV2_N_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TRSYLV2_T_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TRSYLV2_T_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)

    ADD_TEST(SINGLE_TRSYLV2_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV2_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV2_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV2_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)
    ADD_TEST(SINGLE_TRSYLV2_N_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(SINGLE_TRSYLV2_N_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(SINGLE_TRSYLV2_T_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(SINGLE_TRSYLV2_T_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_trsylv2_single --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)

    ENDIF()
ENDFOREACH()

#TGSYLV (SIGN = 1)
FOREACH(SOLVER RANGE 1 24 1)
	IF (SOLVER LESS 13 OR SOLVER GREATER 18 )
        SET(ROWSX ${ROWS})
        SET(COLSX ${COLS})
    ELSE()
        SET(ROWSX ${ROWS_BIG})
        SET(COLSX ${COLS_BIG})
    ENDIF()

    ADD_TEST(DOUBLE_TGSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TGSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TGSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TGSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)

    ADD_TEST(DOUBLE_TGSYLV_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TGSYLV_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TGSYLV_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(DOUBLE_TGSYLV_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)


    # Role change of B and D
    ADD_TEST(DOUBLE_TGSYLV_N_N_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1 --changerole)
	ADD_TEST(DOUBLE_TGSYLV_N_T_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1 --changerole)
	ADD_TEST(DOUBLE_TGSYLV_T_N_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1 --changerole)
	ADD_TEST(DOUBLE_TGSYLV_T_T_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1 --changerole)

    ADD_TEST(DOUBLE_TGSYLV_N_N_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1 --changerole)
	ADD_TEST(DOUBLE_TGSYLV_N_T_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1 --changerole)
	ADD_TEST(DOUBLE_TGSYLV_T_N_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1 --changerole)
	ADD_TEST(DOUBLE_TGSYLV_T_T_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1 --changerole)


    # Sign change
	ADD_TEST(DOUBLE_TGSYLV_N_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TGSYLV_N_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TGSYLV_T_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TGSYLV_T_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)
    ADD_TEST(DOUBLE_TGSYLV_N_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TGSYLV_N_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TGSYLV_T_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(DOUBLE_TGSYLV_T_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)


ENDFOREACH()


FOREACH(SOLVER RANGE 1 24 1)
	IF (SOLVER LESS 13 OR SOLVER GREATER 18 )
        SET(ROWSX ${ROWS})
        SET(COLSX ${COLS})
    ELSE()
        SET(ROWSX ${ROWS_BIG})
        SET(COLSX ${COLS_BIG})
    ENDIF()

	ADD_TEST(SINGLE_TGSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TGSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TGSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TGSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)

    ADD_TEST(SINGLE_TGSYLV_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TGSYLV_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TGSYLV_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=1)
	ADD_TEST(SINGLE_TGSYLV_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=1)


	ADD_TEST(SINGLE_TGSYLV_N_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TGSYLV_N_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TGSYLV_T_N_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TGSYLV_T_T_SIGN_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)

    ADD_TEST(SINGLE_TGSYLV_N_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TGSYLV_N_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TGSYLV_T_N_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign=-1)
	ADD_TEST(SINGLE_TGSYLV_T_T_SIGN_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign=-1)


ENDFOREACH()

#TCGSYLV (SIGN1 = 1, SIGN2 = 1)
FOREACH(SOLVER RANGE 1 24 1)
	IF (SOLVER LESS 13 OR SOLVER GREATER 18 )
        SET(ROWSX ${ROWS})
        SET(COLSX ${COLS})
    ELSE()
        SET(ROWSX ${ROWS_BIG})
        SET(COLSX ${COLS_BIG})
    ENDIF()



    ADD_TEST(DOUBLE_TGCSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGCSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGCSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGCSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(DOUBLE_TGCSYLV_N_N_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --changerole)
	ADD_TEST(DOUBLE_TGCSYLV_N_T_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --changerole)
	ADD_TEST(DOUBLE_TGCSYLV_T_N_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --changerole)
	ADD_TEST(DOUBLE_TGCSYLV_T_T_BD_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --changerole)

    IF( SOLVER LESS 25 )
    ADD_TEST(DOUBLE_TGCSYLV_N_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_N_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_T_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_T_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(DOUBLE_TGCSYLV_N_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(DOUBLE_TGCSYLV_N_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(DOUBLE_TGCSYLV_T_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(DOUBLE_TGCSYLV_T_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(DOUBLE_TGCSYLV_N_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_N_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_T_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_T_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ENDIF()

    #
    # Small
    #
    ADD_TEST(DOUBLE_TGCSYLV_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGCSYLV_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGCSYLV_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGCSYLV_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(DOUBLE_TGCSYLV_N_N_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --changerole)
    ADD_TEST(DOUBLE_TGCSYLV_N_T_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --changerole)
    ADD_TEST(DOUBLE_TGCSYLV_T_N_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --changerole)
    ADD_TEST(DOUBLE_TGCSYLV_T_T_BD_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --changerole)

    IF( SOLVER LESS 25 )
    ADD_TEST(DOUBLE_TGCSYLV_N_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_N_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_T_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_T_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(DOUBLE_TGCSYLV_N_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(DOUBLE_TGCSYLV_N_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(DOUBLE_TGCSYLV_T_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(DOUBLE_TGCSYLV_T_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(DOUBLE_TGCSYLV_N_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_N_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_T_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_T_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ENDIF()

ENDFOREACH()


FOREACH(SOLVER RANGE 1 24 1)
	IF (SOLVER LESS 13 OR SOLVER GREATER 18 )
        SET(ROWSX ${ROWS})
        SET(COLSX ${COLS})
    ELSE()
        SET(ROWSX ${ROWS_BIG})
        SET(COLSX ${COLS_BIG})
    ENDIF()



    ADD_TEST(SINGLE_TGCSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGCSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGCSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGCSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(SINGLE_TGCSYLV_N_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_N_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_T_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_T_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(SINGLE_TGCSYLV_N_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(SINGLE_TGCSYLV_N_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(SINGLE_TGCSYLV_T_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(SINGLE_TGCSYLV_T_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(SINGLE_TGCSYLV_N_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_N_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_T_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_T_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)


    #
    # Small
    #

    ADD_TEST(SINGLE_TGCSYLV_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TGCSYLV_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TGCSYLV_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TGCSYLV_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(SINGLE_TGCSYLV_N_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_N_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_T_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_T_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(SINGLE_TGCSYLV_N_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(SINGLE_TGCSYLV_N_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(SINGLE_TGCSYLV_T_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(SINGLE_TGCSYLV_T_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(SINGLE_TGCSYLV_N_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_N_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_T_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_T_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)


    # TGCSYLV_DUAL
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_S1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_S-1_S1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single  --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
	ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_S-1_S-1_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWSX} --cols=${COLSX} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)


    #SMALL
    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_N_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(SINGLE_TGCSYLV_DUAL_T_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_S1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=1 --sign2=-1)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_S-1_S1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=1)

    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_N_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_N_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)
    ADD_TEST(DOUBLE_TGCSYLV_DUAL_T_T_S-1_S-1_${SOLVER}_SMALL ${CMAKE_CURRENT_BINARY_DIR}/benchmark_tgcsylv_dual   --nmat=${NMAT} --rows=${ROWS_SMALL} --cols=${COLS_SMALL} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T --solver=${SOLVER} --alignoff --sign1=-1 --sign2=-1)




ENDFOREACH()

