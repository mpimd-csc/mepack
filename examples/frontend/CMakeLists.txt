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
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_DGELYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGELYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_SGELYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGELYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gelyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()


# GESTEIN
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_DGESTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_SGESTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gestein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

# GGLYAP
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_DGGLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGGLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_SGGLYAP_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGGLYAP_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gglyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

# GGSTEIN
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_DGGSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGGSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 24)
	ADD_TEST(FRONTEND_SGGSTEIN_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGGSTEIN_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()


#
# test with m!=n for sylv
#
SET(COLS 169)
# GESYLV
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27)
	ADD_TEST(FRONTEND_DGESYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27)
	ADD_TEST(FRONTEND_SGESYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()


# GESYLV2
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
	ADD_TEST(FRONTEND_DGESYLV2_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESYLV2_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESYLV2_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGESYLV2_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
	ADD_TEST(FRONTEND_SGESYLV2_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESYLV2_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESYLV2_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGESYLV2_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()

# GGSYLV
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27)
	ADD_TEST(FRONTEND_DGGSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGGSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGGSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGGSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27)
	ADD_TEST(FRONTEND_SGGSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGGSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGGSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGGSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()
# GCSYLV
FOREACH(SOLVER 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
	ADD_TEST(FRONTEND_DGCSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGCSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGCSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGCSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

	ADD_TEST(FRONTEND_SGCSYLV_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGCSYLV_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGCSYLV_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGCSYLV_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

    ADD_TEST(FRONTEND_DGCSYLV_DUAL_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGCSYLV_DUAL_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGCSYLV_DUAL_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_DGCSYLV_DUAL_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

	ADD_TEST(FRONTEND_SGCSYLV_DUAL_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGCSYLV_DUAL_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGCSYLV_DUAL_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_TEST(FRONTEND_SGCSYLV_DUAL_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/frontend_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)


ENDFOREACH()



