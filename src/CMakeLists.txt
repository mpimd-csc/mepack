INCLUDE(GenerateExportHeader)

SET (DOUBLE_KERNEL_SRC
	    double/kernels/dla_tgsylv_kron_kernel.f90
	    double/kernels/dla_tgsylv_kernel_44_nn.f90 double/kernels/dla_tgsylv_kernel_44_nt.f90 double/kernels/dla_tgsylv_kernel_44_tn.f90 double/kernels/dla_tgsylv_kernel_44_tt.f90
	    double/kernels/dla_trsylv_kernel_44_nn.f90  double/kernels/dla_trsylv_kernel_44_nt.f90 double/kernels/dla_trsylv_kernel_44_tn.f90 double/kernels/dla_trsylv_kernel_44_tt.f90
	    double/kernels/dla_trsylv2_kernel_44_nn.f90  double/kernels/dla_trsylv2_kernel_44_nt.f90 double/kernels/dla_trsylv2_kernel_44_tn.f90 double/kernels/dla_trsylv2_kernel_44_tt.f90
	    double/kernels/dla_trlyap_kernel_44_n.f90   double/kernels/dla_trlyap_kernel_44_t.f90
	    double/kernels/dla_tgcsylv_kernel_44_nn.f90   double/kernels/dla_tgcsylv_kernel_44_nt.f90  double/kernels/dla_tgcsylv_kernel_44_tn.f90 double/kernels/dla_tgcsylv_kernel_44_tt.f90
	    double/kernels/dla_tgcsylv_dual_kernel_44_nn.f90   double/kernels/dla_tgcsylv_dual_kernel_44_nt.f90  double/kernels/dla_tgcsylv_dual_kernel_44_tn.f90
        double/kernels/dla_tgcsylv_dual_kernel_44_tt.f90

)


SET(SINGLE_KERNEL_SRC
	    single/kernels/sla_trsylv_kernel_44_nn.f90  single/kernels/sla_trsylv_kernel_44_nt.f90 single/kernels/sla_trsylv_kernel_44_tn.f90 single/kernels/sla_trsylv_kernel_44_tt.f90
	    single/kernels/sla_trsylv2_kernel_44_nn.f90  single/kernels/sla_trsylv2_kernel_44_nt.f90 single/kernels/sla_trsylv2_kernel_44_tn.f90 single/kernels/sla_trsylv2_kernel_44_tt.f90
	    single/kernels/sla_tgsylv_kernel_44_nn.f90 single/kernels/sla_tgsylv_kernel_44_nt.f90 single/kernels/sla_tgsylv_kernel_44_tn.f90 single/kernels/sla_tgsylv_kernel_44_tt.f90
	    single/kernels/sla_trlyap_kernel_44_n.f90   single/kernels/sla_trlyap_kernel_44_t.f90
	    single/kernels/sla_tgcsylv_kernel_44_nn.f90   single/kernels/sla_tgcsylv_kernel_44_nt.f90  single/kernels/sla_tgcsylv_kernel_44_tn.f90 single/kernels/sla_tgcsylv_kernel_44_tt.f90
        single/kernels/sla_tgcsylv_dual_kernel_44_nn.f90   single/kernels/sla_tgcsylv_dual_kernel_44_nt.f90  single/kernels/sla_tgcsylv_dual_kernel_44_tn.f90
        single/kernels/sla_tgcsylv_dual_kernel_44_tt.f90
        single/kernels/sla_tgsylv_kron_kernel.f90

)
SET (DOUBLE_LEVEL2_SRC
	    double/level2/dla_tgsylv_l2_opt_local_copy_32.f90
        double/level2/dla_tgsylv_l2_opt_local_copy_64.f90
	    double/level2/dla_tgsylv_l2_opt_local_copy_96.f90
	    double/level2/dla_tgsylv_l2_opt_local_copy_128.f90
	    double/level2/dla_tgsylv_l2_opt_local_copy.f90
	    double/level2/dla_tgsylv_l2_reorder.f90
	    double/level2/dla_tgsylv_l2.f90
	    double/level2/dla_tgsylv_gardiner_laub.f90
	    double/level2/dla_tgsylv_l2_rowwise.f90
	    double/level2/dla_tgsylv_l2_colwise.f90
	    double/level2/dla_tglyap_l2.f90
	    double/level2/dla_tglyap_l2_unopt.f90
	    double/level2/dla_tgstein_l2.f90
	    double/level2/dla_trsylv_l2_unopt.f90
	    double/level2/dla_trsylv_l2.f90
	    double/level2/dla_trstein_l2.f90
	    double/level2/dla_trsylv_l2_opt_reorder.f90
	    double/level2/dla_trsylv_l2_opt_local_copy.f90
	    double/level2/dla_trsylv_l2_opt_local_copy_32.f90
	    double/level2/dla_trsylv_l2_opt_local_copy_64.f90
	    double/level2/dla_trsylv_l2_opt_local_copy_96.f90
	    double/level2/dla_trsylv_l2_opt_local_copy_128.f90
	    double/level2/dla_trsylv2_l2_unopt.f90
	    double/level2/dla_trsylv2_l2.f90
	    double/level2/dla_trsylv2_l2_opt_reorder.f90
	    double/level2/dla_trsylv2_l2_opt_local_copy.f90
	    double/level2/dla_trsylv2_l2_opt_local_copy_32.f90
	    double/level2/dla_trsylv2_l2_opt_local_copy_64.f90
	    double/level2/dla_trsylv2_l2_opt_local_copy_96.f90
        double/level2/dla_trsylv2_l2_opt_local_copy_128.f90
	    double/level2/dla_trlyap_l2.f90
	    double/level2/dla_trlyap_l2_opt.f90
	    double/level2/dla_tgcsylv_l2_unopt.f90
	    double/level2/dla_tgcsylv_l2.f90
	    double/level2/dla_tgcsylv_l2_opt_reorder.f90
	    double/level2/dla_tgcsylv_l2_opt_local_copy.f90
	    double/level2/dla_tgcsylv_l2_opt_local_copy_32.f90
	    double/level2/dla_tgcsylv_l2_opt_local_copy_64.f90
	    double/level2/dla_tgcsylv_l2_opt_local_copy_96.f90
	    double/level2/dla_tgcsylv_l2_opt_local_copy_128.f90
        double/level2/dla_tgcsylv_dual_l2.f90
        double/level2/dla_tgcsylv_dual_l2_opt_local_copy.f90
        double/level2/dla_tgcsylv_dual_l2_opt_local_copy_32.f90
        double/level2/dla_tgcsylv_dual_l2_opt_local_copy_64.f90
        double/level2/dla_tgcsylv_dual_l2_opt_local_copy_96.f90
        double/level2/dla_tgcsylv_dual_l2_opt_local_copy_128.f90
        )


SET (SINGLE_LEVEL2_SRC
	    single/level2/sla_tgsylv_l2_opt_local_copy_32.f90
	    single/level2/sla_tgsylv_l2_opt_local_copy_64.f90
	    single/level2/sla_tgsylv_l2_opt_local_copy_96.f90
	    single/level2/sla_tgsylv_l2_opt_local_copy_128.f90
	    single/level2/sla_tgsylv_l2_opt_local_copy.f90
	    single/level2/sla_tgsylv_l2_reorder.f90
	    single/level2/sla_tgsylv_l2.f90
	    single/level2/sla_tgsylv_gardiner_laub.f90
	    single/level2/sla_tgsylv_l2_rowwise.f90
        single/level2/sla_tgsylv_l2_colwise.f90
	    single/level2/sla_tglyap_l2.f90
	    single/level2/sla_tglyap_l2_unopt.f90
	    single/level2/sla_tgstein_l2.f90
	    single/level2/sla_trsylv_l2_unopt.f90
	    single/level2/sla_trsylv_l2.f90
	    single/level2/sla_trstein_l2.f90
	    single/level2/sla_trsylv_l2_opt_reorder.f90
	    single/level2/sla_trsylv_l2_opt_local_copy.f90
	    single/level2/sla_trsylv_l2_opt_local_copy_32.f90
	    single/level2/sla_trsylv_l2_opt_local_copy_64.f90
	    single/level2/sla_trsylv_l2_opt_local_copy_96.f90
	    single/level2/sla_trsylv_l2_opt_local_copy_128.f90
	    single/level2/sla_trsylv2_l2_unopt.f90
	    single/level2/sla_trsylv2_l2.f90
	    single/level2/sla_trsylv2_l2_opt_reorder.f90
	    single/level2/sla_trsylv2_l2_opt_local_copy.f90
	    single/level2/sla_trsylv2_l2_opt_local_copy_32.f90
	    single/level2/sla_trsylv2_l2_opt_local_copy_64.f90
	    single/level2/sla_trsylv2_l2_opt_local_copy_96.f90
	    single/level2/sla_trsylv2_l2_opt_local_copy_128.f90
	    single/level2/sla_trlyap_l2.f90
	    single/level2/sla_trlyap_l2_opt.f90
	    single/level2/sla_tgcsylv_l2_unopt.f90
	    single/level2/sla_tgcsylv_l2.f90
	    single/level2/sla_tgcsylv_l2_opt_reorder.f90
	    single/level2/sla_tgcsylv_l2_opt_local_copy.f90
	    single/level2/sla_tgcsylv_l2_opt_local_copy_32.f90
	    single/level2/sla_tgcsylv_l2_opt_local_copy_64.f90
	    single/level2/sla_tgcsylv_l2_opt_local_copy_96.f90
	    single/level2/sla_tgcsylv_l2_opt_local_copy_128.f90
        single/level2/sla_tgcsylv_dual_l2.f90
        single/level2/sla_tgcsylv_dual_l2_opt_local_copy.f90
        single/level2/sla_tgcsylv_dual_l2_opt_local_copy_32.f90
        single/level2/sla_tgcsylv_dual_l2_opt_local_copy_64.f90
        single/level2/sla_tgcsylv_dual_l2_opt_local_copy_96.f90
        single/level2/sla_tgcsylv_dual_l2_opt_local_copy_128.f90
        )

SET (DOUBLE_LEVEL3_SRC
	    double/level3/dla_tgsylv_l3_opt.f90
	    double/level3/dla_tgsylv_l3_rowwise.f90
	    double/level3/dla_tgsylv_l3_colwise.f90

	    double/level3/dla_tglyap_l3.f90
	    double/level3/dla_tgstein_l3.f90
	    double/level3/dla_trsylv_l3_unopt.f90
	    double/level3/dla_trsylv_l3_opt.f90
	    double/level3/dla_trsylv2_l3_opt.f90
	    double/level3/dla_trsylv2_l3_unopt.f90
	    double/level3/dla_trlyap_l3.f90
	    double/level3/dla_trstein_l3.f90
	    double/level3/dla_tgcsylv_l3_unopt.f90
	    double/level3/dla_tgcsylv_l3_opt.f90
	    double/level3/dla_tgcsylv_dual_l3.f90


        )

SET (SINGLE_LEVEL3_SRC
	    single/level3/sla_tgsylv_l3_opt.f90
	    single/level3/sla_tgsylv_l3_rowwise.f90
        single/level3/sla_tgsylv_l3_colwise.f90
	    single/level3/sla_tglyap_l3.f90
	    single/level3/sla_tgstein_l3.f90
	    single/level3/sla_trsylv_l3_unopt.f90
	    single/level3/sla_trsylv_l3_opt.f90
	    single/level3/sla_trsylv2_l3_opt.f90
	    single/level3/sla_trsylv2_l3_unopt.f90
	    single/level3/sla_trlyap_l3.f90
	    single/level3/sla_trstein_l3.f90
	    single/level3/sla_tgcsylv_l3_unopt.f90
	    single/level3/sla_tgcsylv_l3_opt.f90
	    single/level3/sla_tgcsylv_dual_l3.f90


	)

SET (DOUBLE_RECURSIVE_SRC
	    double/recursive/dla_trlyap_recursive.f90
	    double/recursive/dla_trstein_recursive.f90
	    double/recursive/dla_trsylv2_recursive.f90
	    double/recursive/dla_trsylv_recursive.f90
	    double/recursive/dla_tgsylv_recursive.f90
	    double/recursive/dla_tgcsylv_recursive.f90
	    double/recursive/dla_tglyap_recursive.f90
	    double/recursive/dla_tgstein_recursive.f90
	    double/recursive/dla_tgcsylv_dual_recursive.f90

        )

SET (SINGLE_RECURSIVE_SRC
	    single/recursive/sla_trlyap_recursive.f90
	    single/recursive/sla_trstein_recursive.f90
	    single/recursive/sla_trsylv2_recursive.f90
	    single/recursive/sla_trsylv_recursive.f90
	    single/recursive/sla_tgsylv_recursive.f90
	    single/recursive/sla_tgcsylv_recursive.f90
	    single/recursive/sla_tglyap_recursive.f90
	    single/recursive/sla_tgstein_recursive.f90
	    single/recursive/sla_tgcsylv_dual_recursive.f90


	)

SET (DOUBLE_2STAGE_SRC
	    double/level3/dla_tgstein_l3_2stage.f90
	    double/level3/dla_tglyap_l3_2stage.f90
	    double/level3/dla_trsylv_l3_2stage.f90
	    double/level3/dla_tgsylv_l3_2stage.f90
	    double/level3/dla_trsylv2_l3_2stage.f90
	    double/level3/dla_trlyap_l3_2stage.f90
	    double/level3/dla_trstein_l3_2stage.f90
	    double/level3/dla_tgcsylv_l3_2stage.f90
        double/level3/dla_tgcsylv_dual_l3_2stage.f90


    )
SET (SINGLE_2STAGE_SRC
	    single/level3/sla_tgstein_l3_2stage.f90
	    single/level3/sla_tglyap_l3_2stage.f90
	    single/level3/sla_trsylv_l3_2stage.f90
	    single/level3/sla_tgsylv_l3_2stage.f90
	    single/level3/sla_trsylv2_l3_2stage.f90
	    single/level3/sla_trlyap_l3_2stage.f90
	    single/level3/sla_trstein_l3_2stage.f90
	    single/level3/sla_tgcsylv_l3_2stage.f90
        single/level3/sla_tgcsylv_dual_l3_2stage.f90
        )

SET (DOUBLE_DAG_SRC
	    double/dag/dla_tgsylv_dag.f90
	    double/dag/dla_tgcsylv_dag.f90
	    double/dag/dla_tgstein_dag.f90
	    double/dag/dla_tglyap_dag.f90
	    double/dag/dla_trsylv_dag.f90
	    double/dag/dla_trsylv2_dag.f90
        double/dag/dla_trlyap_dag.f90
	    double/dag/dla_trstein_dag.f90
	    double/dag/dla_tgcsylv_dual_dag.f90

)

SET (SINGLE_DAG_SRC
	    single/dag/sla_tgsylv_dag.f90
	    single/dag/sla_tgcsylv_dag.f90
	    single/dag/sla_tgstein_dag.f90
	    single/dag/sla_tglyap_dag.f90
	    single/dag/sla_trsylv_dag.f90
	    single/dag/sla_trsylv2_dag.f90
	    single/dag/sla_trlyap_dag.f90
	    single/dag/sla_trstein_dag.f90
        single/dag/sla_tgcsylv_dual_dag.f90

        )

SET(DOUBLE_FRONTEND
    double/frontend/dla_ggcsylv.f90
    double/frontend/dla_gelyap.f90
    double/frontend/dla_gestein.f90
    double/frontend/dla_gesylv.f90
    double/frontend/dla_gesylv2.f90
    double/frontend/dla_gglyap.f90
    double/frontend/dla_ggstein.f90
    double/frontend/dla_ggsylv.f90
    double/frontend/dla_ggcsylv_dual.f90


	double/refine/dla_ggsylv_refine.f90
	double/refine/dla_ggcsylv_refine.f90
	double/refine/dla_ggstein_refine.f90
	double/refine/dla_gglyap_refine.f90
	double/refine/dla_gelyap_refine.f90
	double/refine/dla_gestein_refine.f90
    double/refine/dla_gesylv_refine.f90
    double/refine/dla_gesylv2_refine.f90
	double/refine/dla_ggcsylv_dual_refine.f90





)

SET(SINGLE_FRONTEND
    single/frontend/sla_ggcsylv.f90
    single/frontend/sla_gelyap.f90
    single/frontend/sla_gestein.f90
    single/frontend/sla_gesylv.f90
    single/frontend/sla_gesylv2.f90
    single/frontend/sla_gglyap.f90
    single/frontend/sla_ggstein.f90
    single/frontend/sla_ggsylv.f90
    single/frontend/sla_ggcsylv_dual.f90


    single/refine/sla_ggsylv_refine.f90
	single/refine/sla_ggcsylv_refine.f90
	single/refine/sla_ggstein_refine.f90
	single/refine/sla_gglyap_refine.f90
	single/refine/sla_gelyap_refine.f90
	single/refine/sla_gestein_refine.f90
    single/refine/sla_gesylv_refine.f90
    single/refine/sla_gesylv2_refine.f90
	single/refine/sla_ggcsylv_dual_refine.f90
    )

SET(DOUBLE_SRC ${DOUBLE_KERNEL_SRC}  ${DOUBLE_LEVEL2_SRC} ${DOUBLE_LEVEL3_SRC} ${DOUBLE_RECURSIVE_SRC} ${DOUBLE_FRONTEND})
SET(SINGLE_SRC ${SINGLE_KERNEL_SRC}  ${SINGLE_LEVEL2_SRC} ${SINGLE_LEVEL3_SRC} ${SINGLE_RECURSIVE_SRC} ${SINGLE_FRONTEND})
SET(DOUBLE_SRC ${DOUBLE_SRC} ${DOUBLE_2STAGE_SRC} ${DOUBLE_DAG_SRC})
SET(SINGLE_SRC ${SINGLE_SRC} ${SINGLE_2STAGE_SRC} ${SINGLE_DAG_SRC})

SET(MISC_SRC
	   misc/lapack_fix/dgetc2x.f90 misc/lapack_fix/sgetc2x.f90 misc/lapack_fix/cgetc2x.f90 misc/lapack_fix/zgetc2x.f90
	   misc/lapack_fix/dtrmvs.f90 misc/lapack_fix/strmvs.f90
	   misc/dla_qtrmm2.f90 misc/sla_qtrmm2.f90
	   misc/dla_dtrmm_fix.f90 misc/sla_strmm_fix.f90
	   misc/dla_geaxb.f90 misc/sla_geaxb.f90
	   misc/dla_sort_gev.f90 misc/sla_sort_gev.f90
	   misc/dla_sort_ev.f90 misc/sla_sort_ev.f90
       misc/dla_small_solve4.f90 misc/sla_small_solve4.f90
	   misc/dla_small_solve8.f90 misc/sla_small_solve8.f90
	   misc/dla_itranspose.f90 misc/sla_itranspose.f90
	   misc/dla_transform_standard.f90 misc/sla_transform_standard.f90
	   misc/dla_transform_general.f90 misc/sla_transform_general.f90
	   misc/xerror_handler.c

       common/options/machine_double.f90
       common/options/machine_single.f90
       common/options/blocksize_double.f90
       common/options/blocksize_single.f90
       common/options/isolver.f90
       common/options/isolver.c
       common/options/blocksize_double.c
       common/options/blocksize_single.c
       common/options/sortev.f90
       common/options/sortev.c
       common/options/openmp.f90
       common/options/openmp.c
       common/options/frontend_solver.f90
       common/options/frontend_solver.c
       common/options/verbose.f90
       common/options/verbose.c
       common/misc/epsilon.c
       common/misc/memory.c
       common/misc/memory_frontend.c
       common/init.f90
       common/init.c
       common/init_lua.c
       common/version.c
       common/version.F90
       )

SET (C_INTERFACE_SRC
    c/dag/trlyap.c
    c/level3/trlyap.c
    c/level2/trlyap.c
    c/recursive/trlyap.c
    c/dag/trstein.c
    c/level3/trstein.c
    c/level2/trstein.c
    c/recursive/trstein.c
    c/dag/tglyap.c
    c/level3/tglyap.c
    c/level2/tglyap.c
    c/recursive/tglyap.c
    c/dag/tgstein.c
    c/level3/tgstein.c
    c/level2/tgstein.c
    c/recursive/tgstein.c
    c/dag/trsylv.c
    c/level2/trsylv.c
    c/level3/trsylv.c
    c/recursive/trsylv.c
    c/dag/trsylv2.c
    c/level2/trsylv2.c
    c/level3/trsylv2.c
    c/recursive/trsylv2.c
    c/dag/tgsylv.c
    c/level2/tgsylv.c
    c/level3/tgsylv.c
    c/recursive/tgsylv.c
    c/dag/tgcsylv.c
    c/level2/tgcsylv.c
    c/level3/tgcsylv.c
    c/recursive/tgcsylv.c
    c/frontend/gelyap.c
    c/frontend/gestein.c
    c/frontend/gesylv.c
    c/frontend/gesylv2.c
    c/frontend/gglyap.c
    c/frontend/ggstein.c
    c/frontend/ggsylv.c
    c/frontend/cgsylv.c
    c/refine/gelyap.c
    c/refine/gestein.c
    c/refine/gglyap.c
    c/refine/ggstein.c
    c/refine/gesylv.c
    c/refine/gesylv2.c
    c/refine/ggsylv.c
    c/refine/cgsylv.c
    c/level2/tgcsylv_dual.c
    c/level3/tgcsylv_dual.c
    c/recursive/tgcsylv_dual.c
    c/dag/tgcsylv_dual.c
    c/frontend/cgsylv_dual.c
    c/refine/cgsylv_dual.c
)



SET(SRC ${DOUBLE_SRC} ${SINGLE_SRC} ${C_INTERFACE_SRC} ${MISC_SRC})

# Integrate the version number
SET_SOURCE_FILES_PROPERTIES(common/version.c PROPERTIES COMPILE_DEFINITIONS "MEPACK_MAJOR=${PROJECT_VERSION_MAJOR};MEPACK_MINOR=${PROJECT_VERSION_MINOR};MEPACK_PATCH_LEVEL=${PROJECT_VERSION_PATCH}")
SET_SOURCE_FILES_PROPERTIES(common/version.F90 PROPERTIES COMPILE_DEFINITIONS "MEPACK_MAJOR=${PROJECT_VERSION_MAJOR};MEPACK_MINOR=${PROJECT_VERSION_MINOR};MEPACK_PATCH_LEVEL=${PROJECT_VERSION_PATCH}")

# Add the libraries
ADD_LIBRARY(mepack SHARED ${SRC} $<TARGET_OBJECTS:lapack_missing> $<TARGET_OBJECTS:cscutils_lua_src> $<TARGET_OBJECTS:lua_src>)
ADD_LIBRARY(mepack_static STATIC ${SRC} $<TARGET_OBJECTS:lapack_missing> $<TARGET_OBJECTS:cscutils_lua_src> $<TARGET_OBJECTS:lua_src> )

# Generate the export headers
generate_export_header(mepack
    EXPORT_FILE_NAME mepack_exports.h
    INCLUDE_GUARD_NAME MEPACK${MEPACK_SUFFIX}_EXPORTS_H
    BASE_NAME mepack
)
SET_TARGET_PROPERTIES(mepack_static PROPERTIES COMPILE_FLAGS -DMEPACK_STATIC_DEFINE)


TARGET_INCLUDE_DIRECTORIES(
        mepack
        PUBLIC
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/f90mod-shared>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
            $<INSTALL_INTERFACE:include/mepack${MEPACK_SUFFIX}>  # <prefix>/include
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}/shared>
        PRIVATE
            $<TARGET_PROPERTY:cscutils_src,INCLUDE_DIRECTORIES>
     )

TARGET_INCLUDE_DIRECTORIES(
        mepack_static
        PUBLIC
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/f90mod-static>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
            $<INSTALL_INTERFACE:include/mepack${MEPACK_SUFFIX}>  # <prefix>/include
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}/static>
        PRIVATE
            $<TARGET_PROPERTY:cscutils_src,INCLUDE_DIRECTORIES>
     )



MACRO(MEPACK_LINK_LIBRARIES target)
    TARGET_LINK_LIBRARIES(
        ${target}
            PUBLIC
                OpenMP::OpenMP_Fortran
                ${LAPACK_LIBRARIES}
                ${BLAS_LIBRARIES}
        )
ENDMACRO(MEPACK_LINK_LIBRARIES)

MEPACK_LINK_LIBRARIES(mepack)
MEPACK_LINK_LIBRARIES(mepack_static)



SET_TARGET_PROPERTIES(mepack PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/f90mod-shared")
SET_TARGET_PROPERTIES(mepack_static PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/f90mod-static")

SET_TARGET_PROPERTIES(mepack PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR}  )
SET_PROPERTY(TARGET mepack PROPERTY
                           PUBLIC_HEADER
                                "${CMAKE_CURRENT_SOURCE_DIR}/mepack.h"
                                "${CMAKE_CURRENT_BINARY_DIR}/mepack_exports.h"
                     )
SET_PROPERTY(TARGET mepack_static PROPERTY
                                  PUBLIC_HEADER
                                     "${CMAKE_CURRENT_SOURCE_DIR}/mepack.h"
                                     "${CMAKE_CURRENT_BINARY_DIR}/mepack_exports.h"
                     )

SET_TARGET_PROPERTIES(mepack PROPERTIES OUTPUT_NAME mepack${MEPACK_SUFFIX})
SET_TARGET_PROPERTIES(mepack_static PROPERTIES OUTPUT_NAME mepack${MEPACK_SUFFIX})



#
# Default LUA Config
#
INCLUDE(${PROJECT_SOURCE_DIR}/libcscutils/src/lua/LuaInterface.cmake)
IF ( NOT MEPACK_LUA_CONFIG )
    CSC_LUA_ADD_INLINE_CODE(mepack ${CMAKE_CURRENT_SOURCE_DIR}/default_config.lua VARNAME mepack_default_config_code FUNCTIONNAME mepack_default_config)
    CSC_LUA_ADD_INLINE_CODE(mepack_static ${CMAKE_CURRENT_SOURCE_DIR}/default_config.lua VARNAME mepack_default_config_code FUNCTIONNAME mepack_default_config)
ELSE()
    CSC_LUA_ADD_INLINE_CODE(mepack ${MEPACK_LUA_CONFIG} VARNAME mepack_default_config_code FUNCTIONNAME mepack_default_config)
    CSC_LUA_ADD_INLINE_CODE(mepack_static ${MEPACK_LUA_CONFIG} VARNAME mepack_default_config_code FUNCTIONNAME mepack_default_config)
ENDIF()



# Fix missing files in coverage
IF (COVERAGE STREQUAL ON)
    FILE(COPY
        double/level2/dla_tgcsylv_l2_opt_local_copy_xx.f90
        double/level2/dla_tgsylv_l2_opt_local_copy_xx.f90
        double/level2/dla_trsylv_l2_opt_local_copy_xx.f90
        double/level2/dla_trsylv2_l2_opt_local_copy_xx.f90
        double/level2/dla_tgcsylv_dual_l2_opt_local_copy_xx.f90

        single/level2/sla_tgcsylv_l2_opt_local_copy_xx.f90
        single/level2/sla_tgcsylv_dual_l2_opt_local_copy_xx.f90
        single/level2/sla_tgsylv_l2_opt_local_copy_xx.f90
        single/level2/sla_trsylv_l2_opt_local_copy_xx.f90
        single/level2/sla_trsylv2_l2_opt_local_copy_xx.f90
        double/level2/inc_solve8.f90
        DESTINATION
        ${CMAKE_CURRENT_BINARY_DIR}
    )
ENDIF()


ADD_LIBRARY(MEPACK::mepack ALIAS mepack)
ADD_LIBRARY(MEPACK::mepack_static ALIAS mepack_static)



INSTALL(TARGETS mepack_static mepack
    EXPORT MEPACKTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mepack${MEPACK_SUFFIX}

)

INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/f90mod-static/" DESTINATION ${CMAKE_INSTALL_MODULEDIR}/static)
INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/f90mod-shared/" DESTINATION ${CMAKE_INSTALL_MODULEDIR}/shared)



IF( MATLAB STREQUAL ON AND NOT MEXOCT_LINK_STATIC STREQUAL ON)
    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/matlab")
    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/matlab/octave")

    ADD_CUSTOM_COMMAND(TARGET mepack POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_SONAME_FILE:mepack>" "${PROJECT_BINARY_DIR}/matlab/matlab"
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_SONAME_FILE:mepack>" "${PROJECT_BINARY_DIR}/matlab/octave")
ENDIF()
