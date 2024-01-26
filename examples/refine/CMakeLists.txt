SET( SRC
    refine_ggsylv
	refine_ggcsylv
	refine_ggstein
	refine_gglyap
    refine_gelyap
    refine_gestein
    refine_gesylv
    refine_gesylv2
    refine_ggcsylv_dual

    refine_gelyap_single
	refine_gestein_single
    refine_gesylv_single
    refine_gesylv2_single
    refine_gglyap_single
    refine_ggstein_single
    refine_ggsylv_single
    refine_ggcsylv_single
    refine_ggcsylv_dual_single


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


SET(ROWS 257)
SET(COLS 169)
SET(MB 32)
SET(NB 64)
SET(NMAT 1)


FOREACH(SOLVER 1 2 3 4 5)
    ADD_MEPACK_TEST(FRONTEND_DGGSYLV_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGGSYLV_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGGSYLV_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGGSYLV_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGGSYLV_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGGSYLV_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGGSYLV_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGGSYLV_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 )
    ADD_MEPACK_TEST(FRONTEND_DGCSYLV_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGCSYLV_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGCSYLV_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGCSYLV_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGCSYLV_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGCSYLV_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGCSYLV_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGCSYLV_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

    ADD_MEPACK_TEST(FRONTEND_DGCSYLV_DUAL_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGCSYLV_DUAL_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGCSYLV_DUAL_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_DGCSYLV_DUAL_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGCSYLV_DUAL_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGCSYLV_DUAL_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGCSYLV_DUAL_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
	ADD_MEPACK_TEST(FRONTEND_SGCSYLV_DUAL_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggcsylv_dual_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)



ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5 )
    ADD_MEPACK_TEST(FRONTEND_DGGSTEIN_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGGSTEIN_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggstein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGGSTEIN_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGGSTEIN_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_ggstein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)

ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5)
    ADD_MEPACK_TEST(FRONTEND_DGGLYAP_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGGLYAP_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGGLYAP_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGGLYAP_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)

ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5)
    ADD_MEPACK_TEST(FRONTEND_DGELYAP_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGELYAP_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGELYAP_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap_single  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGELYAP_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gelyap_single  --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5)
    ADD_MEPACK_TEST(FRONTEND_DGESTEIN_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gestein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESTEIN_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gestein   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESTEIN_REFINE_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gestein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESTEIN_REFINE_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gestein_single   --nmat=${NMAT} --rows=${ROWS}  --mb=${MB}  --runs=1 --trans=T  --solver=${SOLVER} --alignoff)
ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

ENDFOREACH()

FOREACH(SOLVER 1 2 3 4 5)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV2_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV2_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV2_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_DGESYLV2_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV2_REFINE_N_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV2_REFINE_T_N_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=N  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV2_REFINE_N_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=N --transb=T  --solver=${SOLVER} --alignoff)
    ADD_MEPACK_TEST(FRONTEND_SGESYLV2_REFINE_T_T_${SOLVER} ${CMAKE_CURRENT_BINARY_DIR}/refine_gesylv2_single   --nmat=${NMAT} --rows=${ROWS} --cols=${COLS} --mb=${MB} --nb=${NB} --runs=1 --transa=T --transb=T  --solver=${SOLVER} --alignoff)

ENDFOREACH()




