SET (SRC
	double/check_arguments_trsylv
	double/check_arguments_trsylv2
	double/check_arguments_trlyap
	double/check_arguments_trstein
	double/check_arguments_tglyap
	double/check_arguments_tgstein
	double/check_arguments_tgsylv
	double/check_arguments_tcgsylv

    double/check_arguments_gelyap
	double/check_arguments_gestein
	double/check_arguments_gesylv
	double/check_arguments_gesylv2
	double/check_arguments_gglyap
	double/check_arguments_ggstein
	double/check_arguments_ggsylv
	double/check_arguments_cgsylv

    double/check_arguments_ggsylv_refine
	double/check_arguments_cgsylv_refine
	double/check_arguments_gstein_refine
	double/check_arguments_glyap_refine
    double/check_arguments_gelyap_refine
    double/check_arguments_gestein_refine
    double/check_arguments_gesylv_refine
    double/check_arguments_gesylv2_refine

    single/check_arguments_single_cgsylv
    single/check_arguments_single_ggsylv
    single/check_arguments_single_gesylv
    single/check_arguments_single_gesylv2
    single/check_arguments_single_gelyap
    single/check_arguments_single_gestein
    single/check_arguments_single_gglyap
    single/check_arguments_single_ggstein

    single/check_arguments_single_trsylv
	single/check_arguments_single_trsylv2
	single/check_arguments_single_trlyap
	single/check_arguments_single_trstein
	single/check_arguments_single_tglyap
	single/check_arguments_single_tgstein
	single/check_arguments_single_tgsylv
	single/check_arguments_single_tcgsylv

	single/check_arguments_single_ggsylv_refine
	single/check_arguments_single_cgsylv_refine
	single/check_arguments_single_gstein_refine
	single/check_arguments_single_glyap_refine
    single/check_arguments_single_gelyap_refine
    single/check_arguments_single_gestein_refine
    single/check_arguments_single_gesylv_refine
    single/check_arguments_single_gesylv2_refine
    )

foreach (srcfile ${SRC})
    get_filename_component(SRCTGT ${srcfile} NAME)
    add_executable (${SRCTGT} ${srcfile}.c)
    target_link_libraries (${SRCTGT} mepack benchmark_lib cscutils)
    target_include_directories(${SRCTGT} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach (srcfile)

ADD_TEST(ARGUMENTS_TRSYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_trsylv)
ADD_TEST(ARGUMENTS_TRSYLV2 ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_trsylv2)
ADD_TEST(ARGUMENTS_TGSYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_tgsylv)
ADD_TEST(ARGUMENTS_TCGSYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_tcgsylv)
ADD_TEST(ARGUMENTS_TRLYAP ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_trlyap)
ADD_TEST(ARGUMENTS_TGLYAP ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_tglyap)
ADD_TEST(ARGUMENTS_TRSTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_trstein)
ADD_TEST(ARGUMENTS_TGSTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_tgstein)

ADD_TEST(ARGUMENTS_GELYAP ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gelyap)
ADD_TEST(ARGUMENTS_GESTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gestein)
ADD_TEST(ARGUMENTS_GESYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gesylv)
ADD_TEST(ARGUMENTS_GESYLV2 ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gesylv2)
ADD_TEST(ARGUMENTS_GGLYAP ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gglyap)
ADD_TEST(ARGUMENTS_GGSTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_ggstein)
ADD_TEST(ARGUMENTS_GGSYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_ggsylv)
ADD_TEST(ARGUMENTS_CGSYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_cgsylv)

ADD_TEST(ARGUMENTS_GGSYLV_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_ggsylv_refine)
ADD_TEST(ARGUMENTS_CGSYLV_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_cgsylv_refine)
ADD_TEST(ARGUMENTS_GSTEIN_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gstein_refine)
ADD_TEST(ARGUMENTS_GLYAP_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_glyap_refine)
ADD_TEST(ARGUMENTS_GELYAP_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gelyap_refine)
ADD_TEST(ARGUMENTS_GESTEIN_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gestein_refine)
ADD_TEST(ARGUMENTS_GESYLV_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gesylv_refine)
ADD_TEST(ARGUMENTS_GESYLV2_REFINE ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_gesylv2_refine)

ADD_TEST(ARGUMENTS_SINGLE_TRSYLV  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_trsylv)
ADD_TEST(ARGUMENTS_SINGLE_TRSYLV2 ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_trsylv2)
ADD_TEST(ARGUMENTS_SINGLE_TGSYLV  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_tgsylv)
ADD_TEST(ARGUMENTS_SINGLE_TCGSYLV ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_tcgsylv)
ADD_TEST(ARGUMENTS_SINGLE_TRLYAP  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_trlyap)
ADD_TEST(ARGUMENTS_SINGLE_TGLYAP  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_tglyap)
ADD_TEST(ARGUMENTS_SINGLE_TRSTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_trstein)
ADD_TEST(ARGUMENTS_SINGLE_TGSTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_tgstein)

ADD_TEST(ARGUMENTS_SINGLE_GELYAP  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gelyap)
ADD_TEST(ARGUMENTS_SINGLE_GESTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gestein)
ADD_TEST(ARGUMENTS_SINGLE_GESYLV  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gesylv)
ADD_TEST(ARGUMENTS_SINGLE_GESYLV2 ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gesylv2)
ADD_TEST(ARGUMENTS_SINGLE_GGLYAP  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gglyap)
ADD_TEST(ARGUMENTS_SINGLE_GGSTEIN ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_ggstein)
ADD_TEST(ARGUMENTS_SINGLE_GGSYLV  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_ggsylv)
ADD_TEST(ARGUMENTS_SINGLE_CGSYLV  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_cgsylv)

ADD_TEST(ARGUMENTS_SINGLE_GGSYLV_REFINE   ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_ggsylv_refine)
ADD_TEST(ARGUMENTS_SINGLE_CGSYLV_REFINE   ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_cgsylv_refine)
ADD_TEST(ARGUMENTS_SINGLE_GSTEIN_REFINE   ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gstein_refine)
ADD_TEST(ARGUMENTS_SINGLE_GLYAP_REFINE    ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_glyap_refine)
ADD_TEST(ARGUMENTS_SINGLE_GELYAP_REFINE   ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gelyap_refine)
ADD_TEST(ARGUMENTS_SINGLE_GESTEIN_REFINE  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gestein_refine)
ADD_TEST(ARGUMENTS_SINGLE_GESYLV_REFINE   ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gesylv_refine)
ADD_TEST(ARGUMENTS_SINGLE_GESYLV2_REFINE  ${CMAKE_CURRENT_BINARY_DIR}/check_arguments_single_gesylv2_refine)















