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

#ifndef MEPACK_INTERNAL_H

#define MEPACK_INTERNAL_H


#if __GNUC__ > 7
#define FORTRAN_LEN_T size_t
#else
#define FORTRAN_LEN_T int
#endif

#ifdef INTEGER8
#include <stdint.h>
#define Int int64_t
#endif
#ifndef Int
#define Int int
#endif

/* Initialization */
void FC_GLOBAL_(mepack_initialize,MEPACK_INITIALIZE) (void);

/*
 * Blocksizes
 */
void FC_MODULE_(mepack_options_blocksize_double,trsylv_blocksize_mb_set,MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trsylv_blocksize_nb_set,MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_double,trsylv2_blocksize_mb_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV2_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trsylv2_blocksize_nb_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV2_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_double,tgsylv_blocksize_mb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSYLV_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgsylv_blocksize_nb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSYLV_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_double,tgcsylv_blocksize_mb_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgcsylv_blocksize_nb_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_double,tgcsylv_dual_blocksize_mb_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_DUAL_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgcsylv_dual_blocksize_nb_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_DUAL_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_double,trlyap_blocksize_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRLYAP_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tglyap_blocksize_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGLYAP_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trstein_blocksize_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSTEIN_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgstein_blocksize_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSTEIN_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trlyap_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRLYAP_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tglyap_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGLYAP_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trstein_blocksize_2stage_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSTEIN_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgstein_blocksize_2stage_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSTEIN_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trsylv_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,trsylv2_blocksize_2stage_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV2_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgsylv_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSYLV_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgcsylv_blocksize_2stage_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_double,tgcsylv_dual_blocksize_2stage_set, MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_DUAL_BLOCKSIZE_2STAGE_SET)( Int *MB);

void FC_MODULE_(mepack_options_blocksize_single,trsylv_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trsylv_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_single,trsylv2_blocksize_mb_set, MEPACK_BLOCKSIZE_SINGLE,TRSYLV2_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trsylv2_blocksize_nb_set, MEPACK_BLOCKSIZE_SINGLE,TRSYLV2_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_single,tgsylv_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TGSYLV_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgsylv_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TGSYLV_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_single,tgcsylv_blocksize_mb_set, MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgcsylv_blocksize_nb_set, MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_BLOCKSIZE_NB_SET)( Int *NB);
void FC_MODULE_(mepack_options_blocksize_single,tgcsylv_dual_blocksize_mb_set, MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_DUAL_BLOCKSIZE_MB_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgcsylv_dual_blocksize_nb_set, MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_DUAL_BLOCKSIZE_NB_SET)( Int *NB);


void FC_MODULE_(mepack_options_blocksize_single,trlyap_blocksize_set,  MEPACK_BLOCKSIZE_SINGLE,TRLYAP_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tglyap_blocksize_set,  MEPACK_BLOCKSIZE_SINGLE,TGLYAP_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trstein_blocksize_set, MEPACK_BLOCKSIZE_SINGLE,TRSTEIN_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgstein_blocksize_set, MEPACK_BLOCKSIZE_SINGLE,TGSTEIN_BLOCKSIZE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trlyap_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TRLYAP_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tglyap_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGLYAP_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trstein_blocksize_2stage_set, MEPACK_BLOCKSIZE_SINGLE,TRSTEIN_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgstein_blocksize_2stage_set, MEPACK_BLOCKSIZE_SINGLE,TGSTEIN_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trsylv_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,trsylv2_blocksize_2stage_set, MEPACK_BLOCKSIZE_SINGLE,TRSYLV2_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgsylv_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGSYLV_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgcsylv_blocksize_2stage_set, MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_BLOCKSIZE_2STAGE_SET)( Int *MB);
void FC_MODULE_(mepack_options_blocksize_single,tgcsylv_dual_blocksize_2stage_set, MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_DUAL_BLOCKSIZE_2STAGE_SET)( Int *MB);


/*
 * Frontends
 */
void FC_MODULE_(mepack_options_frontend_solver,trsylv_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRSYLV2_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,trsylv2_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRSYLV22_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,tgsylv_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGSYLV_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,tgcsylv_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGCSYLV_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,tgcsylv_dual_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGCSYLV_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,trlyap_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRLYAP_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,trstein_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRSTEIN_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,tglyap_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGLYAP_FRONTEND_SOLVER_SET)( Int *FSOLVE);
void FC_MODULE_(mepack_options_frontend_solver,tgstein_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGSTEIN_FRONTEND_SOLVER_SET)( Int *FSOLVE);

/*
 * Isolver
 */
void FC_MODULE_(mepack_options_isolver,trsylv_isolver_set,MEPACK_OPTIONS_ISOLVER,TRSYLV2_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,trsylv2_isolver_set,MEPACK_OPTIONS_ISOLVER,TRSYLV22_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,tgsylv_isolver_set,MEPACK_OPTIONS_ISOLVER,TGSYLV_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,tgcsylv_isolver_set,MEPACK_OPTIONS_ISOLVER,TGCSYLV_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,tgcsylv_dual_isolver_set,MEPACK_OPTIONS_ISOLVER,TGCSYLV_ISOLVER_SET)( Int *ISOLVE);

void FC_MODULE_(mepack_options_isolver,trlyap_isolver_set,MEPACK_OPTIONS_ISOLVER,TRLYAP_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,trstein_isolver_set,MEPACK_OPTIONS_ISOLVER,TRSTEIN_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,tglyap_isolver_set,MEPACK_OPTIONS_ISOLVER,TGLYAP_ISOLVER_SET)( Int *ISOLVE);
void FC_MODULE_(mepack_options_isolver,tgstein_isolver_set,MEPACK_OPTIONS_ISOLVER,TGSTEIN_ISOLVER_SET)( Int *ISOLVE);

/*
 * SortEV
 */
Int  FC_MODULE_(mepack_options_sortev,sortev_enabled,MEPACK_OPTIONS_SORTEV,SORTEV_ENABLED)(void);
void FC_MODULE_(mepack_options_sortev,sortev_enable,MEPACK_OPTIONS_SORTEV,SORTEV_ENABLE)(Int* SORT);

/*
 * OpenMP
 */
Int  FC_MODULE_(mepack_options_openmp,openmp_enabled,MEPACK_OPTIONS_OPENMP,OPENMP_ENABLED)(void);
void FC_MODULE_(mepack_options_openmp,openmp_enable,MEPACK_OPTIONS_OPENMP,OPENMP_ENABLE)(Int* openmp);

/*
 * Verbosity
 */
void FC_MODULE_(mepack_options_verbose,verbose_init,MEPACK_OPTIONS_VERBOSE,VERBOSE_INIT)(void);
void FC_MODULE_(mepack_options_verbose,verbose_set,MEPACK_OPTIONS_VERBOSE,VERBOSE_SET)(Int level);
Int  FC_MODULE_(mepack_options_verbose,verbose_get,MEPACK_OPTIONS_VERBOSE,VERBOSE_GET)(void);

/*
 *  TRLYAP (DOUBLE)
 */
void FC_GLOBAL_(dla_trlyap_dag,DLA_TRLYAP_DAG)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trlyap_l3,DLA_TRLYAP_L3)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trlyap_l3_2s,DLA_TRLYAP_L3_2S)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trlyap_l2,DLA_TRLYAP_L2)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T  l1);
void FC_GLOBAL_(dla_trlyap_l2_opt,DLA_TRLYAP_L2_OPT)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T  l1);
void FC_GLOBAL_(dla_trlyap_recursive,DLA_TRLYAP_RECURSIVE)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);

/*
 *  TRLYAP (SINGLE)
 */
void FC_GLOBAL_(sla_trlyap_dag,SLA_TRLYAP_DAG)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trlyap_l3,SLA_TRLYAP_L3)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trlyap_l3_2s,SLA_TRLYAP_L3_2S)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trlyap_l2,SLA_TRLYAP_L2)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trlyap_l2_opt,SLA_TRLYAP_L2_OPT)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trlyap_recursive,SLA_TRLYAP_RECURSIVE)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
/*
 * TRSTEIN  (DOUBLE)
 */
void FC_GLOBAL_(dla_trstein_l3,DLA_TRSTEIN_L3)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trstein_dag,DLA_TRSTEIN_DAG)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trstein_l3_2s,DLA_TRSTEIN_L3_2S)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trstein_recursive,DLA_TRSTEIN_RECURSIVE)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_trstein_l2,DLA_TRSTEIN_L2)(const char *TRANSA, Int *M, double * A, Int *LDA, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);

/*
 * TRSTEIN (SINGLE)
 */
void FC_GLOBAL_(sla_trstein_l3,SLA_TRSTEIN_L3)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trstein_dag,SLA_TRSTEIN_DAG)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trstein_l3_2s,SLA_TRSTEIN_L3_2S)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trstein_recursive,SLA_TRSTEIN_RECURSIVE)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_trstein_l2,SLA_TRSTEIN_L2)(const char *TRANSA, Int *M, float * A, Int *LDA, float *X, Int *LDX, float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);

/*
 * TRSYLV (DOUBLE)
 */
void FC_GLOBAL_(dla_trsylv_l3,DLA_TRSYLV_L3)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l3_unopt,DLA_TRSYLV_L3_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_trsylv_l3_2s,DLA_TRSYLV_L3_2S)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_dag,DLA_TRSYLV_DAG)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_trsylv_l2,DLA_TRSYLV_L2)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_reorder,DLA_TRSYLV_L2_REORDER)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_unopt,DLA_TRSYLV_L2_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_local_copy,DLA_TRSYLV_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_recursive,DLA_TRSYLV_RECURSIVE)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_local_copy_32,DLA_TRSYLV_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_local_copy_64,DLA_TRSYLV_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_local_copy_96,DLA_TRSYLV_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv_l2_local_copy_128,DLA_TRSYLV_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * TRSYLV (SINGLE)
 */
void FC_GLOBAL_(sla_trsylv_l3,SLA_TRSYLV_L3)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l3_unopt,SLA_TRSYLV_L3_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(sla_trsylv_l3_2s,SLA_TRSYLV_L3_2S)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_dag,SLA_TRSYLV_DAG)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(sla_trsylv_l2,SLA_TRSYLV_L2)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_reorder,SLA_TRSYLV_L2_REORDER)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_unopt,SLA_TRSYLV_L2_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_local_copy,SLA_TRSYLV_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_recursive,SLA_TRSYLV_RECURSIVE)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_local_copy_32,SLA_TRSYLV_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_local_copy_64,SLA_TRSYLV_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_local_copy_96,SLA_TRSYLV_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv_l2_local_copy_128,SLA_TRSYLV_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
/*
 * TRSYLV2 (DOUBLE)
 */

void FC_GLOBAL_(dla_trsylv2_l3,DLA_TRSYLV2_L3)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l3_unopt,DLA_TRSYLV2_L3_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_dag,DLA_TRSYLV2_DAG)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l3_2s,DLA_TRSYLV2_L3_2S)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_trsylv2_l2,DLA_TRSYLV2_L2)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_reorder,DLA_TRSYLV2_L2_REORDER)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_unopt,DLA_TRSYLV2_L2_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_local_copy,DLA_TRSYLV2_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_recursive,DLA_TRSYLV2_RECURSIVE)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_local_copy_32,DLA_TRSYLV2_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_local_copy_64,DLA_TRSYLV2_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_local_copy_96,DLA_TRSYLV2_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_trsylv2_l2_local_copy_128,DLA_TRSYLV2_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * TRSYLV2 (SINGLE)
 */
void FC_GLOBAL_(sla_trsylv2_l3,SLA_TRSYLV2_L3)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l3_unopt,SLA_TRSYLV2_L3_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_dag,SLA_TRSYLV2_DAG)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l3_2s,SLA_TRSYLV2_L3_2S)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(sla_trsylv2_l2,SLA_TRSYLV2_L2)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_reorder,SLA_TRSYLV2_L2_REORDER)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_unopt,SLA_TRSYLV2_L2_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_local_copy,SLA_TRSYLV2_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_recursive,SLA_TRSYLV2_RECURSIVE)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_local_copy_32,SLA_TRSYLV2_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_local_copy_64,SLA_TRSYLV2_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_local_copy_96,SLA_TRSYLV2_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_trsylv2_l2_local_copy_128,SLA_TRSYLV2_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
/*
 *  TGLYAP(DOUBLE)
 */
void FC_GLOBAL_(dla_tglyap_l3,DLA_TGLYAP_L3)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1 );
void FC_GLOBAL_(dla_tglyap_l2,DLA_TGLYAP_L2)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tglyap_l2_unopt,DLA_TGLYAP_L2_UNOPT)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tglyap_dag,DLA_TGLYAP_DAG)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tglyap_l3_2s,DLA_TGLYAP_L3_2S)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tglyap_recursive,DLA_TGLYAP_RECURSIVE)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);


/*
 * TGLYAP(SINGLE)
 */

void FC_GLOBAL_(sla_tglyap_l3,SLA_TGLYAP_L3)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tglyap_l2,SLA_TGLYAP_L2)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tglyap_l2_unopt,SLA_TGLYAP_L2_UNOPT)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tglyap_dag,SLA_TGLYAP_DAG)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tglyap_l3_2s,SLA_TGLYAP_L3_2S)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tglyap_recursive,SLA_TGLYAP_RECURSIVE)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);

/*
 * TGSTEIN(DOUBLE)
 */
void FC_GLOBAL_(dla_tgstein_l3,DLA_TGSTEIN_L3)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tgstein_dag,DLA_TGSTEIN_DAG)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tgstein_l3_2s,DLA_TGSTEIN_L3_2S)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tgstein_recursive,DLA_TGSTEIN_RECURSIVE)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(dla_tgstein_l2,DLA_TGSTEIN_L2)(const char *TRANSA, Int *M, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);

/*
 * TGSTEIN(SINGLE)
 */
void FC_GLOBAL_(sla_tgstein_l3,SLA_TGSTEIN_L3)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tgstein_dag,SLA_TGSTEIN_DAG)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tgstein_l3_2s,SLA_TGSTEIN_L3_2S)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tgstein_recursive,SLA_TGSTEIN_RECURSIVE)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);
void FC_GLOBAL_(sla_tgstein_l2,SLA_TGSTEIN_L2)(const char *TRANSA, Int *M, float * A, Int *LDA,float * B, Int * LDB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1);




/*
 * Solvers for TGSYLV (Double)
 */
void FC_GLOBAL_(dla_tgsylv_l3,DLA_TGSYLV_L3)(const char *TRANSA, const char*TRANSB, double *sgn,
        Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l3_rowwise,DLA_TGSYLV_L3_ROWWISE)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l3_colwise,DLA_TGSYLV_L3_COLWISE)(const char *TRANSA, const char*TRANSB, double *sgn,
        Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_tgsylv_l3_2s,DLA_TGSYLV_L3_2S)(const char *TRANSA, const char*TRANSB, double *sgn,
        Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_dag,DLA_TGSYLV_DAG)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_tgsylv_l2_local_copy_32,DLA_TGSYLV_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_local_copy_64,DLA_TGSYLV_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_local_copy_96,DLA_TGSYLV_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_local_copy_128,DLA_TGSYLV_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_local_copy,DLA_TGSYLV_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_reorder,DLA_TGSYLV_L2_REORDER)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_tgsylv_l2,DLA_TGSYLV_L2)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_rowwise,DLA_TGSYLV_L2_ROWWISE)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_l2_colwise,DLA_TGSYLV_L2_COLWISE)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_tgsylv_l2_unopt,DLA_TGSYLV_L2_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(dla_tgsylv_recursive,DLA_TGSYLV_RECURSIVE)(const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgsylv_gardiner_laub,DLA_TGSYLV_GARDINER_LAUB)(const char *TRANSA, const char*TRANSB, double *sgn, Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


/*
 * Solvers for TGSYLV (Single)
 */
void FC_GLOBAL_(sla_tgsylv_l3_2s,SLA_TGSYLV_L3_2S)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_dag,SLA_TGSYLV_DAG)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(sla_tgsylv_l3,SLA_TGSYLV_L3)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l3_unopt,SLA_TGSYLV_L3_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

void FC_GLOBAL_(sla_tgsylv_l2,SLA_TGSYLV_L2)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_unopt,SLA_TGSYLV_L2_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_local_copy,SLA_TGSYLV_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_local_copy_32,SLA_TGSYLV_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_local_copy_64,SLA_TGSYLV_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_local_copy_96,SLA_TGSYLV_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_local_copy_128,SLA_TGSYLV_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_reorder,SLA_TGSYLV_L2_REORDER)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_recursive,SLA_TGSYLV_RECURSIVE)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_gardiner_laub,SLA_TGSYLV_GARDINER_LAUB)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float *A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l3_rowwise,SLA_TGSYLV_L3_ROWWISE)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l3_colwise,SLA_TGSYLV_L3_COLWISE)(const char *TRANSA, const char*TRANSB, float *sgn,
        Int *M, Int *N, float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_rowwise,SLA_TGSYLV_L2_ROWWISE)(const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgsylv_l2_colwise,SLA_TGSYLV_L2_COLWISE)(const char *TRANSA, const char*TRANSB, float *sgn,
        Int *M, Int *N, float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * TGCSYLV (DOUBLE)
 */
void FC_GLOBAL_(dla_tgcsylv_l2,DLA_TGCSYLV_L2)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_unopt,DLA_TGCSYLV_L2_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_reorder,DLA_TGCSYLV_L2_REORDER)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_local_copy,DLA_TGCSYLV_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_local_copy_32,DLA_TGCSYLV_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_local_copy_64,DLA_TGCSYLV_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_local_copy_96,DLA_TGCSYLV_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l2_local_copy_128,DLA_TGCSYLV_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l3_unopt,DLA_TGCSYLV_L3_UNOPT)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l3,DLA_TGCSYLV_L3)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_recursive,DLA_TGCSYLV_RECURSIVE)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dag,DLA_TGCSYLV_DAG)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_l3_2s,DLA_TGCSYLV_L3_2S)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * TGCSYLV (SINGLE)
 */
void FC_GLOBAL_(sla_tgcsylv_l2,SLA_TGCSYLV_L2)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_unopt,SLA_TGCSYLV_L2_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_reorder,SLA_TGCSYLV_L2_REORDER)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_local_copy,SLA_TGCSYLV_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_local_copy_32,SLA_TGCSYLV_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_local_copy_64,SLA_TGCSYLV_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_local_copy_96,SLA_TGCSYLV_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l2_local_copy_128,SLA_TGCSYLV_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l3_unopt,SLA_TGCSYLV_L3_UNOPT)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l3,SLA_TGCSYLV_L3)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_recursive,SLA_TGCSYLV_RECURSIVE)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dag,SLA_TGCSYLV_DAG)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_l3_2s,SLA_TGCSYLV_L3_2S)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * TGCSYLV_DUAL (DOUBLE)
 */
void FC_GLOBAL_(dla_tgcsylv_dual_l2,DLA_TGCSYLV_DUAL_L2)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_32,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_64,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_96,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_128,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF , double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l3,DLA_TGCSYLV_DUAL_L3)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_recursive,DLA_TGCSYLV_DUAL_RECURSIVE)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_dag,DLA_TGCSYLV_DUAL_DAG)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(dla_tgcsylv_dual_l3_2s,DLA_TGCSYLV_DUAL_L3_2S)(const char *TRANSA, const char*TRANSB, double *sgn1, double *sgn2,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD,
        double *E, Int *LDE, double *F, Int *LDF ,  double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * TGCSYLV (SINGLE)
 */
void FC_GLOBAL_(sla_tgcsylv_dual_l2,SLA_TGCSYLV_DUAL_L2)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_32,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_64,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_96,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_128,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF , float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l3,SLA_TGCSYLV_DUAL_L3)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_recursive,SLA_TGCSYLV_DUAL_RECURSIVE)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_dag,SLA_TGCSYLV_DUAL_DAG)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_tgcsylv_dual_l3_2s,SLA_TGCSYLV_DUAL_L3_2S)(const char *TRANSA, const char*TRANSB, float *sgn1, float *sgn2,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *C, Int *LDC, float *D, Int *LDD,
        float *E, Int *LDE, float *F, Int *LDF ,  float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);




/*
 * GELYAP
 */
void FC_GLOBAL_(dla_gelyap,DLA_GELYAP)(const char * FACT, const char *TRANSA, Int *M, double * A, Int *LDA, double * Q, Int *LDQ, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *LDWORK, Int *INFO,  FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_gelyap,SLA_GELYAP)(const char * FACT, const char *TRANSA, Int *M, float * A, Int *LDA, float * Q, Int *LDQ, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
/*
 * GSTEIN
 */
void FC_GLOBAL_(dla_gestein,DLA_GELYAP)(const char * FACT, const char *TRANSA, Int *M, double * A, Int *LDA, double * Q, Int *LDQ, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_gestein,SLA_GELYAP)(const char * FACT, const char *TRANSA, Int *M, float * A, Int *LDA, float * Q, Int *LDQ, float *X,
        Int *LDX, float * SCALE, float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
/*
 * GESYLV
 */
void FC_GLOBAL_(dla_gesylv,DLA_GESYLV)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *QA, Int *LDQA, double *QB, Int *LDQB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);
void FC_GLOBAL_(sla_gesylv,SLA_GESYLV)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *QA, Int *LDQA, float *QB, Int *LDQB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);

/*
 * GESYLV2
 */
void FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N,
        double * A, Int *LDA,double * B, Int * LDB, double *QA, Int *LDQA, double *QB, Int *LDQB, double *X, Int *LDX,
        double * SCALE, double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);
void FC_GLOBAL_(sla_gesylv2,SLA_GESYLV2)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float *sgn,  Int *M, Int *N,
        float * A, Int *LDA,float * B, Int * LDB, float *QA, Int *LDQA, float *QB, Int *LDQB, float *X, Int *LDX,
        float * SCALE, float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);

/*
 * CGSYLV
 */
void FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        double *sgn1, double *sgn2,  Int *M, Int *N, double * A, Int *LDA, double * B, Int * LDB,
        double *C, Int *LDC, double *D, Int *LDD, double *QA, Int *LDQA, double *ZA, Int *LDZA,
        double *QB, Int *LDQB, double *ZB, Int *LDZB, double *E, Int *LDE, double *F, Int *LDF, double * SCALE,
        double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);
void FC_GLOBAL_(sla_ggcsylv,SLA_GGCSYLV)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        float *sgn1, float *sgn2,  Int *M, Int *N, float * A, Int *LDA, float * B, Int * LDB,
        float *C, Int *LDC, float *D, Int *LDD, float *QA, Int *LDQA, float *ZA, Int *LDZA,
        float *QB, Int *LDQB, float *ZB, Int *LDZB, float *E, Int *LDE, float *F, Int *LDF, float * SCALE,
        float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);

/*
 * CGSYLV DUAL
 */
void FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        double *sgn1, double *sgn2,  Int *M, Int *N, double * A, Int *LDA, double * B, Int * LDB,
        double *C, Int *LDC, double *D, Int *LDD, double *QA, Int *LDQA, double *ZA, Int *LDZA,
        double *QB, Int *LDQB, double *ZB, Int *LDZB, double *E, Int *LDE, double *F, Int *LDF, double * SCALE,
        double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);
void FC_GLOBAL_(sla_ggcsylv_dual,SLA_GGCSYLV_DUAL)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        float *sgn1, float *sgn2,  Int *M, Int *N, float * A, Int *LDA, float * B, Int * LDB,
        float *C, Int *LDC, float *D, Int *LDD, float *QA, Int *LDQA, float *ZA, Int *LDZA,
        float *QB, Int *LDQB, float *ZB, Int *LDZB, float *E, Int *LDE, float *F, Int *LDF, float * SCALE,
        float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);


/*
 * GGSYLV
 */

void FC_GLOBAL_(dla_ggsylv,DLA_GGSYLV)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        double *sgn,  Int *M, Int *N, double * A, Int *LDA, double * B, Int * LDB,
        double *C, Int *LDC, double *D, Int *LDD, double *QA, Int *LDQA, double *ZA, Int *LDZA,
        double *QB, Int *LDQB, double *ZB, Int *LDZB, double *X, Int *LDX, double * SCALE,
        double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);
void FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        float *sgn,  Int *M, Int *N, float * A, Int *LDA, float * B, Int * LDB,
        float *C, Int *LDC, float *D, Int *LDD, float *QA, Int *LDQA, float *ZA, Int *LDZA,
        float *QB, Int *LDQB, float *ZB, Int *LDZB, float *X, Int *LDX, float * SCALE,
        float *WORK, Int *LDWORK, Int *INFO,  FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);

/*
 * GGLYAP
 */
void FC_GLOBAL_(dla_gglyap,DLA_GGLYAP)(const char * FACT, const char *TRANSA, Int *M, double * A, Int *LDA,
        double *B, Int *LDB, double * Q, Int *LDQ, double *Z, Int *LDZ, double *X, Int *LDX, double * SCALE,
        double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_gglyap,SLA_GGLYAP)(const char * FACT, const char *TRANSA, Int *M, float * A, Int *LDA,
        float *B, Int *LDB, float * Q, Int *LDQ, float *Z, Int *LDZ, float *X, Int *LDX, float * SCALE,
        float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * GGSTEIN
 */
void FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)(const char * FACT, const char *TRANSA, Int *M, double * A, Int *LDA,
        double *B, Int *LDB, double * Q, Int *LDQ, double *Z, Int *LDZ, double *X, Int *LDX, double * SCALE,
        double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_ggstein,SLA_GGSTEIN)(const char * FACT, const char *TRANSA, Int *M, float * A, Int *LDA,
        float *B, Int *LDB, float * Q, Int *LDQ, float *Z, Int *LDZ, float *X, Int *LDX, float * SCALE,
        float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);

/*
 * Refine
 */

/*
 * GELYAP
 */
void FC_GLOBAL_(dla_gelyap_refine, DLA_GELYAP_REFINE) (const char * TRANS, const char * GUESS, Int * M , double * A, Int * LDA,
        double *X, Int *LDX, double *Y, Int *LDY,  double *AS, Int * LDAS, double * Q, Int *LDQ,
        Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_gelyap_refine, SLA_GELYAP_REFINE) (const char * TRANS, const char * GUESS, Int * M , float * A, Int * LDA,
        float *X, Int *LDX, float *Y, Int *ldy, float *AS, Int * LDAS, float * Q, Int *LDQ,
        Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
/*
 * GESTEIN
 */
void FC_GLOBAL_(dla_gestein_refine, DLA_GESTEIN_REFINE) (const char * TRANS, const char * GUESS, Int * M , double * A, Int * LDA,
        double *X, Int *LDX, double *Y, Int *LDY, double *AS, Int * LDAS, double * Q, Int *LDQ,
        Int * MAXIT, double *TAU,double *CONVLOG,  double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_gestein_refine, SLA_GESTEIN_REFINE) (const char * TRANS, const char *GUESS, Int * M , float * A, Int * LDA,
        float *X, Int *LDX, float *Y, Int *LDY, float *AS, Int * LDAS, float * Q, Int *LDQ,
        Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


/*
 * GESYLV
 */
void FC_GLOBAL_(dla_gesylv_refine, DLA_GESYLV_REFINE) (const char * TRANSA, const char * TRANSB, const char *GUESS, double *SGN, Int * M , Int* N,  double * A, Int * LDA,
        double *B, Int *LDB, double *X, Int *LDX, double *Y, Int *LDY,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * Q, Int *LDQ,
        double *U, Int *LDU, Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);
void FC_GLOBAL_(sla_gesylv_refine, SLA_GESYLV_REFINE) (const char * TRANSA, const char * TRANSB,const char *GUESS,  float *SGN, Int * M , Int* N,  float * A, Int * LDA,
        float *B, Int *LDB, float *X, Int *LDX, float *Y, Int *LDY,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * Q, Int *LDQ,
        float *U, Int *LDU, Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);

/*
 * GESYLV2
 */
void FC_GLOBAL_(dla_gesylv2_refine, DLA_GESYLV2_REFINE) (const char * TRANSA, const char * TRANSB, const char *GUESS, double *SGN, Int * M , Int* N,  double * A, Int * LDA,
        double *B, Int *LDB, double *X, Int *LDX, double *Y, Int *LDY,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * Q, Int *LDQ,
        double *U, Int *LDU, Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);
void FC_GLOBAL_(sla_gesylv2_refine, SLA_GESYLV2_REFINE) (const char * TRANSA, const char * TRANSB,const char *GUESS,  float *SGN, Int * M , Int* N,  float * A, Int * LDA,
        float *B, Int *LDB, float *X, Int *LDX, float *Y, Int *LDY,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * Q, Int *LDQ,
        float *U, Int *LDU, Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2 , FORTRAN_LEN_T l3);

/*
 * GGLYAP
 */
void FC_GLOBAL_(dla_gglyap_refine, DLA_GGLYAP_REFINE) (const char * TRANS, const char *GUESS, Int * M , double * A, Int * LDA,
        double *B, Int *LDB, double *X, Int *LDX, double *Y, Int *LDY,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * Q, Int *LDQ, double * Z, Int * LDZ,
        Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_gglyap_refine, SLA_GGLYAP_REFINE) (const char * TRANS, const char *GUESS, Int * M , float * A, Int * LDA,
        float *B, Int *LDB, float *X, Int *LDX, float *Y, Int *LDY,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * Q, Int *LDQ, float * Z, Int * LDZ,
        Int * MAXIT, float *TAU,float *CONVLOG,  float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


/*
 * GGSTEIN
 */
void FC_GLOBAL_(dla_ggstein_refine, DLA_GGSTEIN_REFINE) (const char * TRANS, const char *GUESS, Int * M , double * A, Int * LDA,
        double *B, Int *LDB, double *X, Int *LDX, double *Y, Int *LDY,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * Q, Int *LDQ, double * Z, Int * LDZ,
        Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);
void FC_GLOBAL_(sla_ggstein_refine, SLA_GGSTEIN_REFINE) (const char * TRANS, const char *GUESS, Int * M , float * A, Int * LDA,
        float *B, Int *LDB, float *X, Int *LDX, float *Y, Int *LDY,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * Q, Int *LDQ, float * Z, Int * LDZ,
        Int * MAXIT, float *TAU,float *CONVLOG,  float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


/*
 * GGSYLV
 */
void FC_GLOBAL_(dla_ggsylv_refine, DLA_GGSYLV_REFINE) (const char * TRANSA, const char * TRANSB, const char *GUESS, double *SGN, Int * M , Int* N,  double * A, Int * LDA,
        double *B, Int *LDB, double *C, Int * LDC, double * D, Int * LDD, double *X, Int *LDX, double *Y, Int *LDY,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * CS, Int * LDCS, double * DS, Int * LDDS, double * Q, Int *LDQ, double * Z, Int * LDZ,
        double *U, Int *LDU, double * V, Int * LDV,  Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);
void FC_GLOBAL_(sla_ggsylv_refine, SLA_GGSYLV_REFINE) (const char * TRANSA, const char * TRANSB, const char *GUESS, float *SGN, Int * M , Int* N,  float * A, Int * LDA,
        float *B, Int *LDB, float *C, Int * LDC, float * D, Int * LDD, float *X, Int *LDX, float *Y, Int *LDY,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * CS, Int * LDCS, float * DS, Int * LDDS, float * Q, Int *LDQ, float * Z, Int * LDZ,
        float *U, Int *LDU, float * V, Int * LDV,  Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);

/*
 * GCSYLV
 */
void FC_GLOBAL_(dla_ggcsylv_refine, DLA_GGCSYLV_REFINE) (const char * TRANSA, const char * TRANSB,const char *GUESS,  double *SGN1, double *SGN2, Int * M , Int* N,  double * A, Int * LDA,
        double *B, Int *LDB, double *C, Int * LDC, double * D, Int * LDD, double *R, Int *LDR, double *L, Int *LDL, double *E, Int *LDE, double *F, Int *LDF,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * CS, Int * LDCS, double * DS, Int * LDDS, double * Q, Int *LDQ, double * Z, Int * LDZ,
        double *U, Int *LDU, double * V, Int * LDV,  Int * MAXIT, double *TAU,double *CONVLOG,  double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);
void FC_GLOBAL_(sla_ggcsylv_refine, SLA_GGCSYLV_REFINE) (const char * TRANSA, const char * TRANSB,const char *GUESS,  float *SGN1, float *SGN2, Int * M , Int* N,  float * A, Int * LDA,
        float *B, Int *LDB, float *C, Int * LDC, float * D, Int * LDD, float *R, Int *LDR, float *L, Int *LDL, float *E, Int *LDE, float *F, Int *LDF,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * CS, Int * LDCS, float * DS, Int * LDDS, float * Q, Int *LDQ, float * Z, Int * LDZ,
        float *U, Int *LDU, float * V, Int * LDV,  Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);

/*
 * GCSYLV_DUAL
 */
void FC_GLOBAL_(dla_ggcsylv_dual_refine, DLA_GGCSYLV_DUAL_REFINE) (const char * TRANSA, const char * TRANSB,const char *GUESS,  double *SGN1, double *SGN2, Int * M , Int* N,  double * A, Int * LDA,
        double *B, Int *LDB, double *C, Int * LDC, double * D, Int * LDD, double *R, Int *LDR, double *L, Int *LDL, double *E, Int *LDE, double *F, Int *LDF,
        double *AS, Int * LDAS, double * BS, Int *LDBS, double * CS, Int * LDCS, double * DS, Int * LDDS, double * Q, Int *LDQ, double * Z, Int * LDZ,
        double *U, Int *LDU, double * V, Int * LDV,  Int * MAXIT, double *TAU,double *CONVLOG,  double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);
void FC_GLOBAL_(sla_ggcsylv_dual_refine, SLA_GGCSYLV_DUAL_REFINE) (const char * TRANSA, const char * TRANSB, const char *GUESS, float *SGN1, float *SGN2, Int * M , Int* N,  float * A, Int * LDA,
        float *B, Int *LDB, float *C, Int * LDC, float * D, Int * LDD, float *R, Int *LDR, float *L, Int *LDL, float *E, Int *LDE, float *F, Int *LDF,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * CS, Int * LDCS, float * DS, Int * LDDS, float * Q, Int *LDQ, float * Z, Int * LDZ,
        float *U, Int *LDU, float * V, Int * LDV,  Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);



#endif /* end of include guard: MEPACK_INTERNAL_H */

