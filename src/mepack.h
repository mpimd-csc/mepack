/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2017-2022
 */

#ifndef MEPACK_H

#define MEPACK_H

#include <stdlib.h>
#include <stdint.h>
#include "mepack_exports.h"

#ifdef __cplusplus
extern "C" {
#endif

    /* Machine number information */
    MEPACK_EXPORT double mepack_double_epsilon(void);
    MEPACK_EXPORT float mepack_single_epsilon(void);


    MEPACK_EXPORT void mepack_version(int *major, int *minor, int *patchlevel);
    MEPACK_EXPORT void mepack_init(void);

    /* Memory query */
    MEPACK_EXPORT ssize_t mepack_memory(const char *func, int M, int N);
    MEPACK_EXPORT ssize_t mepack_memory_frontend(const char *func, const char * factA, const char *factB, int M, int N);

    /**
     * @addtogroup frontend_solver
     * @{
     */

    /**
     * @brief Use the default setting for the level-2 inner solver.
     *
     * By default MEPACK uses the LOCAL_COPY solver with alignment. This default
     * may change if more improved solvers are developed.
     */
    #define MEPACK_ISOLVER_DEFAULT 1

    /**
     * @brief Use the level-2 solver with aligned local copies.
     */
    #define MEPACK_ISOLVER_LOCAL_COPY_ALIGN 1

    /**
     * @brief Use the level-2 solver with local copies without alignment.
     */
    #define MEPACK_ISOLVER_LOCAL_COPY 2

    /**
     * @brief Use the level-2 solver with reordered operations.
     */
    #define MEPACK_ISOLVER_REORDER 3

    /**
     * @brief Use the level-2 solver with standard BLAS operations.
     */
    #define MEPACK_ISOLVER_LEVEL2  4

    /**
     * @brief Use the naively implemented level-2 solver.
     */
    #define MEPACK_ISOLVER_LEVEL2_NAIV 5

    /**
     * @brief Use the recursive blocking solvers instead of level-2 solvers.
     */
    #define MEPACK_ISOLVER_RECURSIVE 6

    /**
     * @}
     */


    MEPACK_EXPORT void mepack_trsylv_isolver_set(int IS);
    MEPACK_EXPORT void mepack_trsylv2_isolver_set(int IS);
    MEPACK_EXPORT void mepack_tgsylv_isolver_set(int IS);
    MEPACK_EXPORT void mepack_tgcsylv_isolver_set(int IS);
    MEPACK_EXPORT void mepack_tgcsylv_dual_isolver_set(int IS);
    MEPACK_EXPORT void mepack_trlyap_isolver_set(int IS);
    MEPACK_EXPORT void mepack_tglyap_isolver_set(int IS);
    MEPACK_EXPORT void mepack_trstein_isolver_set(int IS);
    MEPACK_EXPORT void mepack_tgstein_isolver_set(int IS);

    /* Functions to set the frontend solver.  */

    /**
     * @addtogroup frontend_solver
     * @{
     */

    /**
     * @brief Use the level-3 triangular solver in the frontend.
     */
    #define MEPACK_FRONTEND_SOLVER_LEVEL3 1

    /**
     * @brief Use the level-2 triangular solver in the frontend.
     */
    #define MEPACK_FRONTEND_SOLVER_LEVEL2 2

    /**
     * @brief Use the OpenMP 4 Directed Acyclic Graph solver in the frontend.
     */
    #define MEPACK_FRONTEND_SOLVER_DAG 3

    /**
     * @brief Use the Two-Stage Level-3/OpenMP 4 solver in the frontend.
     */
    #define MEPACK_FRONTEND_SOLVER_2STAGE 4

    /**
     * @brief Use the Recursive Blocking solver in the frontend.
     */
    #define MEPACK_FRONTEND_SOLVER_RECURSIVE 5

    /**
     * @brief Use the Gardiner-Laub solver in the frontend.
     *
     * The Gardiner-Laub solver is only available for the generalized Sylvester equation.
     */
    #define MEPACK_FRONTEND_SOLVER_GARDINER_LAUB 6

    /**
     * @brief Use the LAPACK Sylvester solver in the frontend.
     *
     * The LAPACK-Sylvester solver (DTRSYL) is only available for the
     * standard Sylvester equation.
     */
    #define MEPACK_FRONTEND_SOLVER_LAPACK 7
    /**
     * @}
     */

    MEPACK_EXPORT void mepack_trsylv_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_trsylv2_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_tgsylv_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_tgcsylv_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_tgcsylv_dual_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_trlyap_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_tglyap_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_trstein_frontend_solver_set(int FS);
    MEPACK_EXPORT void mepack_tgstein_frontend_solver_set(int FS);


    /*  Functions for setting the blocksize  */
    MEPACK_EXPORT void mepack_double_trsylv_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trsylv_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_double_trsylv2_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trsylv2_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_double_tgsylv_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgsylv_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_double_tgcsylv_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgcsylv_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_blocksize_nb_set ( int nb ) ;

    MEPACK_EXPORT void mepack_double_trlyap_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tglyap_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trstein_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgstein_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trlyap_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tglyap_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trstein_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgstein_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trsylv_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_trsylv2_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgsylv_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgcsylv_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_blocksize_2stage_set ( int mb ) ;


    MEPACK_EXPORT void mepack_single_trsylv_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trsylv_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_single_trsylv2_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trsylv2_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_single_tgsylv_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgsylv_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_single_tgcsylv_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgcsylv_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_blocksize_mb_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_blocksize_nb_set ( int nb ) ;
    MEPACK_EXPORT void mepack_single_trlyap_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tglyap_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trstein_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgstein_blocksize_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trlyap_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tglyap_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trstein_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgstein_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trsylv_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_trsylv2_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgsylv_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgcsylv_blocksize_2stage_set ( int mb ) ;
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_blocksize_2stage_set ( int mb ) ;


    /* Eigenvalue sorting */
    MEPACK_EXPORT int mepack_sortev_enabled(void);
    MEPACK_EXPORT void mepack_sortev_enable(int sort);

    /* OpenMP */
    MEPACK_EXPORT int mepack_openmp_enabled(void);
    MEPACK_EXPORT void mepack_openmp_enable(int openmp);

    /* Verbosity */
    MEPACK_EXPORT void mepack_verbose_init(void);
    MEPACK_EXPORT void mepack_verbose_set(int level);
    MEPACK_EXPORT int mepack_verbose_get(void);

    /*
     * Triangular solvers
     */

    /* TRLYAP */
    MEPACK_EXPORT void mepack_single_trlyap_dag(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_dag(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trlyap_level3(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_level3(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trlyap_level3_2stage(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_level3_2stage(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trlyap_level3_2stage(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_level3_2stage(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trlyap_level2(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_level2(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trlyap_level2_opt(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_level2_opt(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trlyap_recursive(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trlyap_recursive(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    /* TRSTEIN */
    MEPACK_EXPORT void mepack_single_trstein_dag(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trstein_dag(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trstein_level3(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trstein_level3(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trstein_level2(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trstein_level2(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trstein_level3_2stage(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trstein_level3_2stage(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    MEPACK_EXPORT void mepack_single_trstein_recursive(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_trstein_recursive(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO);

    /* TGLYAP */
    MEPACK_EXPORT void mepack_double_tglyap_dag(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tglyap_dag(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tglyap_level3(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tglyap_level3(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tglyap_level3_2stage(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tglyap_level3_2stage(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tglyap_level2(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tglyap_level2(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tglyap_level2_unopt(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tglyap_level2_unopt(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tglyap_recursive(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tglyap_recursive(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    /* TGSTEIN */
    MEPACK_EXPORT void mepack_double_tgstein_dag(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgstein_dag(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tgstein_level3(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgstein_level3(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tgstein_level3_2stage(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgstein_level3_2stage(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tgstein_level2(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgstein_level2(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_tgstein_recursive(const char *TRANS, int M, double * A, int LDA, double * B, int LDB, double *X, int LDX, double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgstein_recursive(const char *TRANS, int M, float * A, int LDA, float * B, int LDB, float *X, int LDX, float * SCALE, float *WORK, int *INFO);

    /* TRSYLV */
    MEPACK_EXPORT void mepack_double_trsylv_dag(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_dag(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level3(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                      double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                      double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level3(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                    float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                    float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level3_unopt(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level3_unopt(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level3_2stage(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level3_2stage(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_recursive(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_recursive(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_unopt(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_unopt(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_reorder(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_reorder(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_local_copy(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_local_copy(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_local_copy_32(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_local_copy_32(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_local_copy_64(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_local_copy_64(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_local_copy_96(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_local_copy_96(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv_level2_local_copy_128(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv_level2_local_copy_128(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);
    /* TRSYLV2 */
    MEPACK_EXPORT void mepack_double_trsylv2_dag(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_dag(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level3(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                      double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                      double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level3(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                    float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                    float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level3_unopt(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level3_unopt(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level3_2stage(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level3_2stage(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_recursive(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_recursive(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_unopt(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_unopt(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_reorder(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_reorder(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_local_copy(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_local_copy(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_local_copy_32(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_local_copy_32(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_local_copy_64(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_local_copy_64(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_local_copy_96(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_local_copy_96(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    MEPACK_EXPORT void mepack_double_trsylv2_level2_local_copy_128(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_trsylv2_level2_local_copy_128(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO);

    /* TGSYLV */
    MEPACK_EXPORT void mepack_double_tgsylv_dag(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_dag(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);


    MEPACK_EXPORT void mepack_double_tgsylv_level2_local_copy_32(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_local_copy_32(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_local_copy_64(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_local_copy_64(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_local_copy_96(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_local_copy_96(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_local_copy_128(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_local_copy_128(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_local_copy(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_local_copy(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_reorder(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_reorder(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_rowwise(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_rowwise(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2_colwise(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2_colwise(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level2(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level2(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_gardiner_laub(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_gardiner_laub(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level3_rowwise(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level3_rowwise(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level3_colwise(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level3_colwise(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level3(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level3(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_level3_2stage(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_level3_2stage(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgsylv_recursive(const char *TRANSA, const char*TRANSB, double SGN,
                                     int M, int N, double * A, int LDA, double * B, int LDB,
                                     double *C, int LDC, double *D, int LDD, double *X, int LDX,
                                     double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgsylv_recursive(const char *TRANSA, const char*TRANSB, float SGN,
                                     int M, int N, float * A, int LDA, float * B, int LDB,
                                     float *C, int LDC, float *D, int LDD, float *X, int LDX,
                                     float * SCALE, float *WORK, int *INFO);

    /* TGCSYLV */

    MEPACK_EXPORT void mepack_double_tgcsylv_level3(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level3_unopt(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level3_2stage(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dag(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_unopt(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_reorder(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_local_copy(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_local_copy_32(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_local_copy_64(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_local_copy_96(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_level2_local_copy_128(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_recursive(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level3(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level3_unopt(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level3_2stage(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dag(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_unopt(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_reorder(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_local_copy(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_local_copy_32(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_local_copy_64(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_local_copy_96(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_level2_local_copy_128(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_recursive(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);

    /* TGCSYLV DUAL */

    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level3(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level3_2stage(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_dag(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level2(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level2_local_copy(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level2_local_copy_32(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level2_local_copy_64(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level2_local_copy_96(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_level2_local_copy_128(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_double_tgcsylv_dual_recursive(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
            double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
            double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level3(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level3_2stage(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_dag(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level2(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level2_local_copy(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level2_local_copy_32(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level2_local_copy_64(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level2_local_copy_96(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_level2_local_copy_128(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);
    MEPACK_EXPORT void mepack_single_tgcsylv_dual_recursive(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
            float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
            float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO);


    /* GELYAP */
    MEPACK_EXPORT void mepack_double_gelyap(const char * FACT, const char *TRANSA, int M, double * A,  int LDA, double * Q, int LDQ, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gelyap(const char * FACT, const char *TRANSA, int M, float * A,  int LDA, float * Q, int LDQ, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int *INFO);

    /* GESTEIN */
    MEPACK_EXPORT void mepack_double_gestein(const char * FACT, const char *TRANSA, int M, double * A,  int LDA, double * Q, int LDQ, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gestein(const char * FACT, const char *TRANSA, int M, float * A,  int LDA, float * Q, int LDQ, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int *INFO);

    /* GESYLV */
    MEPACK_EXPORT void mepack_double_gesylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *QA, int LDQA, double *QB, int LDQB, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO);
    MEPACK_EXPORT void mepack_single_gesylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *QA, int LDQA, float *QB, int LDQB, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO);

    /* GESYLV2 */
    MEPACK_EXPORT void mepack_double_gesylv2(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *QA, int LDQA, double *QB, int LDQB, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO);
    MEPACK_EXPORT void mepack_single_gesylv2(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *QA, int LDQA, float *QB, int LDQB, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO);

    /* GGLYAP */

    MEPACK_EXPORT void mepack_double_gglyap(const char * FACT, const char *TRANS, int M, double * A,  int LDA, double *B, int LDB,
        double * Q, int LDQ, double *Z, int LDZ, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gglyap(const char * FACT, const char *TRANS, int M, float * A,  int LDA, float *B, int LDB,
        float * Q, int LDQ, float *Z, int LDZ, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int *INFO);

    /* GGSTEIN */

    MEPACK_EXPORT void mepack_double_ggstein(const char * FACT, const char *TRANS, int M, double * A,  int LDA, double *B, int LDB,
        double * Q, int LDQ, double *Z, int LDZ, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_ggstein(const char * FACT, const char *TRANS, int M, float * A,  int LDA, float *B, int LDB,
        float * Q, int LDQ, float *Z, int LDZ, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int *INFO);


    /* GGSYLV */
    MEPACK_EXPORT void mepack_double_ggsylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
        double *QA, int LDQA, double *ZA, int LDZA, double *QB, int LDQB, double *ZB, int LDZB,
        double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO);
    MEPACK_EXPORT void mepack_single_ggsylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
        float *QA, int LDQA, float *ZA, int LDZA, float *QB, int LDQB, float *ZB, int LDZB,
        float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO);

    /* CGSYLV */
    MEPACK_EXPORT void mepack_double_ggcsylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
        double *QA, int LDQA, double *ZA, int LDZA, double *QB, int LDQB, double *ZB, int LDZB,
        double *E, int LDE, double *F, int LDF,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO);
    MEPACK_EXPORT void mepack_single_ggcsylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
        float *QA, int LDQA, float *ZA, int LDZA, float *QB, int LDQB, float *ZB, int LDZB,
        float *E, int LDE, float *F, int LDF,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO);

    /* CGSYLV DUAL */
    MEPACK_EXPORT void mepack_double_ggcsylv_dual(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
        double *QA, int LDQA, double *ZA, int LDZA, double *QB, int LDQB, double *ZB, int LDZB,
        double *E, int LDE, double *F, int LDF,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO);
    MEPACK_EXPORT void mepack_single_ggcsylv_dual(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
        float *QA, int LDQA, float *ZA, int LDZA, float *QB, int LDQB, float *ZB, int LDZB,
        float *E, int LDE, float *F, int LDF,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO);


    /* REFINE GELYAP */
    MEPACK_EXPORT void mepack_double_gelyap_refine(const char * TRANS, const char *GUESS, int M , double * A, int LDA,
        double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double * Q, int LDQ,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gelyap_refine(const char * TRANS, const char *GUESS, int M , float * A, int LDA,
        float *X, int LDX, float *Y, int LDY, float *AS, int LDAS, float * Q, int LDQ,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO);

    /* REFINE GESTEIN */
    MEPACK_EXPORT void mepack_double_gestein_refine(const char * TRANS, const char *GUESS, int M , double * A, int LDA,
        double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double * Q, int LDQ,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gestein_refine(const char * TRANS, const char *GUESS, int M , float * A, int LDA,
        float *X, int LDX, float *Y, int LDY,  float *AS, int LDAS, float * Q, int LDQ,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO);

    /* REFINE GGLYAP */
    MEPACK_EXPORT void mepack_double_gglyap_refine(const char * TRANS, const char *GUESS, int M , double * A, int LDA,
        double *B, int LDB, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double * Q, int LDQ, double *Z, int LDZ,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gglyap_refine(const char * TRANS, const char *GUESS, int M , float * A, int LDA,
        float *B, int LDB, float *X, int LDX, float *Y, int LDY,  float *AS, int LDAS, float *BS, int LDBS,
        float * Q, int LDQ, float *Z, int LDZ,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO);

    /* REFINE GGSTEIN */
    MEPACK_EXPORT void mepack_double_ggstein_refine(const char * TRANS,  const char *GUESS,int M , double * A, int LDA,
        double *B, int LDB, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double * Q, int LDQ, double *Z, int LDZ,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_ggstein_refine(const char * TRANS, const char *GUESS, int M , float * A, int LDA,
        float *B, int LDB, float *X, int LDX, float *Y, int LDY,  float *AS, int LDAS, float *BS, int LDBS,
        float * Q, int LDQ, float *Z, int LDZ,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO);

    /* REFINE GESYLV */
    MEPACK_EXPORT void mepack_double_gesylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN, int M, int N , double * A, int LDA,
        double *B, int LDB, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double * Q, int LDQ, double *U, int LDU,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gesylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN, int M, int N , float * A, int LDA,
        float *B, int LDB, float *X, int LDX, float *Y, int LDY,  float *AS, int LDAS, float *BS, int LDBS,
        float * Q, int LDQ, float *U, int LDU,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO);


    /* REFINE GESYLV2 */
    MEPACK_EXPORT void mepack_double_gesylv2_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN, int M, int N , double * A, int LDA,
        double *B, int LDB, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double * Q, int LDQ, double *U, int LDU,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_gesylv2_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN, int M, int N , float * A, int LDA,
        float *B, int LDB, float *X, int LDX, float *Y, int LDY,  float *AS, int LDAS, float *BS, int LDBS,
        float * Q, int LDQ, float *U, int LDU,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO);

    /* REFINE GGSYLV */
    MEPACK_EXPORT void mepack_double_ggsylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN, int M, int N , double * A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD, double *X, int LDX, double *Y, int LDY,
        double *AS, int LDAS, double *BS, int LDBS, double *CS, int LDCS, double *DS, int LDDS,
        double * Q, int LDQ, double *Z, int LDZ, double *U, int LDU, double *V, int LDV,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_ggsylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN, int M, int N , float * A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD, float *X, int LDX, float *Y, int LDY,
        float *AS, int LDAS, float *BS, int LDBS, float *CS, int LDCS, float *DS, int LDDS,
        float * Q, int LDQ, float *Z, int LDZ, float *U, int LDU, float *V, int LDV,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LDWORK, int *INFO);


    /* Refine CGSYLV */
    MEPACK_EXPORT void mepack_double_ggcsylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN1, double SGN2, int M, int N , double * A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD, double *R, int ldr, double *L, int ldl, double *E, int LDE, double *F, int LDF,
        double *AS, int LDAS, double *BS, int LDBS, double *CS, int LDCS, double *DS, int LDDS,
        double * Q, int LDQ, double *Z, int LDZ, double *U, int LDU, double *V, int LDV,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_ggcsylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN1, float SGN2, int M, int N , float * A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD, float *R, int ldr, float *L, int ldl, float *E, int LDE, float *F, int LDF,
        float *AS, int LDAS, float *BS, int LDBS, float *CS, int LDCS, float *DS, int LDDS,
        float * Q, int LDQ, float *Z, int LDZ, float *U, int LDU, float *V, int LDV,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LDWORK, int *INFO);

    /* Refine CGSYLV DUAL */
    MEPACK_EXPORT void mepack_double_ggcsylv_dual_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN1, double SGN2, int M, int N , double * A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD,double *R, int ldr, double *L, int ldl,  double *E, int LDE, double *F, int LDF,
        double *AS, int LDAS, double *BS, int LDBS, double *CS, int LDCS, double *DS, int LDDS,
        double * Q, int LDQ, double *Z, int LDZ, double *U, int LDU, double *V, int LDV,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LDWORK, int *INFO);
    MEPACK_EXPORT void mepack_single_ggcsylv_dual_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN1, float SGN2, int M, int N , float * A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD, float *R, int ldr, float *L, int ldl, float *E, int LDE, float *F, int LDF,
        float *AS, int LDAS, float *BS, int LDBS, float *CS, int LDCS, float *DS, int LDDS,
        float * Q, int LDQ, float *Z, int LDZ, float *U, int LDU, float *V, int LDV,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LDWORK, int *INFO);



#ifdef __cplusplus
};
#endif


#endif /* end of include guard: MEPACK_H */


