SET( SRC
	frontend_gelyap
	frontend_gelyap_single
	frontend_gestein
	frontend_gestein_single
	frontend_gesylv
	frontend_gesylv_single
	frontend_gesylv2
	frontend_gesylv2_single
	frontend_gglyap
	frontend_gglyap_single
	frontend_ggstein
	frontend_ggstein_single
	frontend_ggsylv
	frontend_ggsylv_single
	frontend_ggcsylv
	frontend_ggcsylv_single
	frontend_ggcsylv_dual
    frontend_ggcsylv_dual_single


)

foreach (srcfile ${SRC})
	add_executable (${srcfile} ${srcfile}.c)
    target_link_libraries (
        ${srcfile}
            PRIVATE
                OpenMP::OpenMP_C
                mepack
                benchmark_lib
                ${ADDITIONAL_LIBRARIES}
                ${LAPACK_LIBRARIES}
                ${BLAS_LIBRARIES}
                cscutils)
endforeach (srcfile)


#
# Frontend Tests
#

SET(ROWS 257)
SET(COLS 257)
SET(MB 32)
SET(NB 64)
SET(NMAT 1)

#
# Test with M = N for LYAP and stein
#
# GELYAP
FOREACH(TRANS "N" "T")
    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_DGELYAP_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap --hess --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_DGELYAP_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
    ENDFOREACH()

    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_SGELYAP_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_SGELYAP_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap_single  --hess --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
    ENDFOREACH()


    # GESTEIN
    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_DGESTEIN_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_DGESTEIN_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein --hess  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
    ENDFOREACH()

    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_SGESTEIN_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_SGESTEIN_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein_single  --hess --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
    ENDFOREACH()

    # GGLYAP
    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_DGGLYAP_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_DGGLYAP_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap --hess  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
    ENDFOREACH()

    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_SGGLYAP_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_SGGLYAP_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap_single  --hess  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)

    ENDFOREACH()

    # GGSTEIN
    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_DGGSTEIN_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_DGGSTEIN_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein --hess  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)

    ENDFOREACH()

    FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
        ADD_MEPACK_TEST(FRONTEND_SGGSTEIN_${TRANS}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
        ADD_MEPACK_TEST(FRONTEND_SGGSTEIN_${TRANS}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein_single --hess  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=${TRANS}  --solver=${SOLVER} --alignoff)
    ENDFOREACH()
ENDFOREACH()

#
# test with m!=n for sylv
#
FOREACH(TRANSA "N" "T")
    FOREACH(TRANSB "N" "T")
        SET(COLS 169)
        # GESYLV
        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 )
            ADD_MEPACK_TEST(FRONTEND_DGESYLV_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_DGESYLV_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv   --nmat=${NMAT} --hess --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)

        ENDFOREACH()

        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 )
            ADD_MEPACK_TEST(FRONTEND_SGESYLV_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_SGESYLV_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv_single   --nmat=${NMAT} --hess --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)


        ENDFOREACH()


        # GESYLV2
        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
            ADD_MEPACK_TEST(FRONTEND_DGESYLV2_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_DGESYLV2_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2   --hess --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)


        ENDFOREACH()

        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
            ADD_MEPACK_TEST(FRONTEND_SGESYLV2_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_SGESYLV2_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2_single --hess  --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)


        ENDFOREACH()

        # GGSYLV
        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27)
            ADD_MEPACK_TEST(FRONTEND_DGGSYLV_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_DGGSYLV_${TRANSA}_${TRANSA}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv  --hess --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)


        ENDFOREACH()

        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27)
            ADD_MEPACK_TEST(FRONTEND_SGGSYLV_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_SGGSYLV_${TRANSA}_${TRANSB}__HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv_single --hess  --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)


        ENDFOREACH()
        # GCSYLV
        FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
            ADD_MEPACK_TEST(FRONTEND_DGCSYLV_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_DGCSYLV_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv  --hess  --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)

            ADD_MEPACK_TEST(FRONTEND_SGCSYLV_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_SGCSYLV_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_single  --hess --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)

            ADD_MEPACK_TEST(FRONTEND_DGCSYLV_DUAL_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_DGCSYLV_DUAL_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual --hess  --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)

            ADD_MEPACK_TEST(FRONTEND_SGCSYLV_DUAL_${TRANSA}_${TRANSB}_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)
            ADD_MEPACK_TEST(FRONTEND_SGCSYLV_DUAL_${TRANSA}_${TRANSB}_HESS_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual_single --hess  --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=${TRANSA} --transb=${TRANSB}  --solver=${SOLVER} --alignoff)


        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()



