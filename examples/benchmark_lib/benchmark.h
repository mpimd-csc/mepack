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
 * Copyright (C) Martin Koehler, 2017-2023
 */

#ifdef __cplusplus
extern "C" {
#endif
#include <ctype.h>
#include <math.h>
#include "FCMangle.h"


#ifdef INTEGER8
#include <stdint.h>
#define Int int64_t
#endif
#ifndef Int
#define Int int
#endif

#include "mepack.h"

#define TOL_SINGLE 2000

#define PATH_LEN 32*1024
    extern char * hdfstore_path;

#define MAX(A,B) (((A)>(B))?(A):(B))

#if __GNUC__ > 7
    typedef size_t fortran_charlen_t;
#else
    typedef int fortran_charlen_t;
#endif


#ifdef MEPACK_STATIC_DEFINE
#  define MEPACK_EXPORT
#  define MEPACK_NO_EXPORT
#else
#  ifndef MEPACK_EXPORT
#    ifdef mepack_EXPORTS
        /* We are building this library */
#      define MEPACK_EXPORT __declspec(dllexport)
#    else
#warning import
        /* We are using this library */
#      define MEPACK_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef MEPACK_NO_EXPORT
#    define MEPACK_NO_EXPORT 
#  endif
#endif

    void benchmark_init(void);
    void benchmark_exit(void);

    double get_wtime(void);
    double get_ctime(void);

    void benchmark_print_matrix_double(Int M, Int N, double *X, Int LDX);
    void benchmark_print_matrix_float(Int M, Int N, float *X, Int LDX);
    Int gcd(Int a, Int b);


    void solver_name(int is);
    const char * solver_name_str(int is);



    /*-----------------------------------------------------------------------------
     *  Right Hand Sides
     *-----------------------------------------------------------------------------*/
    void benchmark_rhs_sylv_double(char *TRANSA, char *TRANSB, double sgn,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY);
    void benchmark_rhs_sylv2_double(char *TRANSA, char *TRANSB, double sgn,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY);
    void benchmark_rhs_gsylv_double(char *TRANSA, char *TRANSB, double sgn,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double * C, Int LDC, double *D, Int LDD,
            double *X, Int LDX, double *Y, Int LDY) ;
    void benchmark_rhs_glyap_double(char *TRANSA, Int M, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY);
    void benchmark_rhs_stein_double(char *TRANSA, Int M, double * A, Int LDA, double *X, Int LDX, double *Y, Int LDY);
    void benchmark_rhs_lyap_double(char *TRANSA, Int M, double * A, Int LDA, double *X, Int LDX, double *Y, Int LDY);
    void benchmark_rhs_ggcsylv_double(char *TRANSA, char *TRANSB, double sgn1, double sgn2,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double * C, Int LDC, double *D, Int LDD,
            double *R, Int LDR, double *L, Int LDL, double *Y1, Int LDY1, double *Y2, Int LDY2);
     void benchmark_rhs_ggcsylv_dual_double(char *TRANSA, char *TRANSB, double sgn1, double sgn2,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double * C, Int LDC, double *D, Int LDD,
            double *R, Int LDR, double *L, Int LDL, double *Y1, Int LDY1, double *Y2, Int LDY2);

    void benchmark_rhs_gstein_double(char *TRANSA, Int M, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY);

    double benchmark_check_X_double(Int M, Int N, double *X1, Int LDX1, double *X2, Int LDX2);

    void benchmark_rhs_sylv_float(char *TRANSA, char *TRANSB, float sgn,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY);
    void benchmark_rhs_sylv2_float(char *TRANSA, char *TRANSB, float sgn,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY);
    void benchmark_rhs_gsylv_float(char *TRANSA, char *TRANSB, float sgn,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float * C, Int LDC, float *D, Int LDD,
            float *X, Int LDX, float *Y, Int LDY) ;
    void benchmark_rhs_glyap_float(char *TRANSA, Int M, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY);
    void benchmark_rhs_stein_float(char *TRANSA, Int M, float * A, Int LDA, float *X, Int LDX, float *Y, Int LDY);
    void benchmark_rhs_lyap_float(char *TRANSA, Int M, float * A, Int LDA, float *X, Int LDX, float *Y, Int LDY);
    void benchmark_rhs_ggcsylv_float(char *TRANSA, char *TRANSB, float sgn1, float sgn2,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float * C, Int LDC, float *D, Int LDD,
            float *R, Int LDR, float *L, Int LDL, float *Y1, Int LDY1, float *Y2, Int LDY2);
    void benchmark_rhs_gstein_float(char *TRANSA, Int M, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY);
    void benchmark_rhs_ggcsylv_dual_float(char *TRANSA, char *TRANSB, float sgn1, float sgn2,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float * C, Int LDC, float *D, Int LDD,
            float *R, Int LDR, float *L, Int LDL, float *Y1, Int LDY1, float *Y2, Int LDY2);


    float benchmark_check_X_float(Int M, Int N, float *X1, Int LDX1, float *X2, Int LDX2);

    /*-----------------------------------------------------------------------------
     *  EVP
     *-----------------------------------------------------------------------------*/
    void benchmark_random_evp_double(Int M, Int *iseed, double *Ar, double *Q, double *A, Int sort);
    void benchmark_random_evp_float(Int M, Int *iseed, float *Ar, float *Q, float *A, Int sort);
    void benchmark_random_gevp_double(Int M, Int *iseed, double *Ar, double *Cr, double *Q, double *Z, double *Ao, double *Co,  Int sort);
    void benchmark_random_gevp_float(Int M, Int *iseed, float *Ar, float *Cr, float *Q, float *Z, float *Ao, float *Co,  Int sort);


    /*-----------------------------------------------------------------------------
     *  Recsy Helper
     *-----------------------------------------------------------------------------*/
    int benchmark_sylv_setup_to_type(const char *TRANSA, const char*TRANSB, double sign);
    int benchmark_gsylv_setup_to_type(const char *TRANSA, const char*TRANSB, double sign);
    int benchmark_glyap_setup_to_type(const char *TRANSA);
    int benchmark_gstein_setup_to_type(const char *TRANSA);
    int benchmark_ggcsylv_setup_to_type(const char *TRANSA, const char*TRANSB, double sign);

    /*-----------------------------------------------------------------------------
     *  HESS
     *-----------------------------------------------------------------------------*/
    void benchmark_ghess_double(Int N, double *A, double *Q, double *B, double *Z, double *Work, Int ldwork);
    void benchmark_ghess_single(Int N, float *A, float *Q, float *B, float *Z, float *Work, Int ldwork);

    /*-----------------------------------------------------------------------------
     *  LAPACK and MEPACK stuff
     *-----------------------------------------------------------------------------*/
    void FC_GLOBAL(dgees,DGEES)(char *JOBVSL, char * SORT, void * SELECTG, Int *N, double *A, Int *LDA, Int *SDIM, double *ALPHAR, double *ALPHAI,
            double *VSL, Int *LDVSL, double *WORK, Int * LDWORK, Int *BWORK, Int *INFO);
    void FC_GLOBAL(dgges,DGGES)(char *JOBVSL, char*JOBVSRL, char * SORT, void * SELECTG, Int *N, double *A, Int *LDA, double *C, Int *ldc, Int *SDIM, double *ALPHAR, double *ALPHAI, double *betar,
            double *VSL, Int *LDVSL, double *VSR, Int *LDVR,  double *WORK, Int * LDWORK, Int *BWORK, Int *INFO);


    extern void FC_GLOBAL(dlaset,DLASET)(char *UL, Int *M, Int *N, double *alpha, double *beta, double *X, Int *lda, fortran_charlen_t l1);
    extern void FC_GLOBAL(dlacpy,DLACPY)(char *uplo,Int *M, Int *N, double *A, Int *LDA, double *B, Int *LDB, fortran_charlen_t l1);

    extern void FC_GLOBAL(dlarnv,DLARNV)(Int *i, Int *S, Int *N, double * X);
    extern double FC_GLOBAL(dlange,DLANGE)(char *NORM, Int *N, Int *M, double *A, Int *LDA, double *WORK, fortran_charlen_t l1);
    extern void FC_GLOBAL(dgemm,DGEMM)(char *ta, char *tb, Int *M, Int *N, Int *K, double *alpha, double *A, Int *LDA, double *B, Int *LDB, double *beta, double *C, Int *LDC, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(dgehrd, DGEHRD)(Int *N, Int *ILO, Int *IHI, double *A, Int *LDA, double *TAU, double *WORK, Int *LDWORK, Int *Info);
    extern void FC_GLOBAL(dggbal, DGGBAL)(char* JOB, Int *N, double *A, Int *LDA, double *B, Int *LDB, Int *ILO, Int *IHI, double *LSCALE, double *RSCALE, double *WORK, Int *INFO, fortran_charlen_t l1);
    extern void FC_GLOBAL(dgghrd, DGGHRD)(char *COMPQ, char* COMPZ, Int *N, Int *ILO, Int *IHI, double *A, Int *LDA,
		    double *B, Int *LDB, double *Q, Int *LDQ, double *Z, Int *LDZ, Int *INFO, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(dgeqrf, DGEQRF)(Int *M, Int *N, double *A, Int *LDA, double *TAU, double *WORK, Int *LWORK, Int *INFO);
    extern void FC_GLOBAL(dormqr, DORMQR)(char *SIDE, char *TRANS, Int *M, Int *N, Int *K, double *A, Int *LDA, double *TAU, double *C, Int *LDC, double *WORK, Int *LWORK, Int *INFO, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(dorgqr, DORGQR)(Int *M, Int *N, Int *K, double *A, Int *LDA, double *TAU, double *WORK, Int *LWORK, Int *INFO);

    extern float FC_GLOBAL(slange,SLANGE)(char *NORM, Int *N, Int *M, float *A, Int *LDA, float *WORK, fortran_charlen_t l1);
    extern void FC_GLOBAL(sgemm,SGEMM)(char *ta, char *tb, Int *M, Int *N, Int *K, float *alpha, float *A, Int *LDA, float *B, Int *LDB, float *beta, float *C, Int *LDC, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(slacpy,SLACPY)(char *uplo,Int *M, Int *N, float *A, Int *LDA, float *B, Int *LDB, fortran_charlen_t l1);
    void FC_GLOBAL(sgees,SGEES)(char *JOBVSL, char * SORT, void * SELECTG, Int *N, float *A, Int *LDA, Int *SDIM, float *ALPHAR, float *ALPHAI,
            float *VSL, Int *LDVSL, float *WORK, Int * LDWORK, Int *BWORK, Int *INFO);
    extern void FC_GLOBAL(slarnv,SLARNV)(Int *i, Int *S, Int *N, float * X);
    extern void FC_GLOBAL(slaset,SLASET)(char *UL, Int *M, Int *N, float *alpha, float *beta, float *X, Int *lda, fortran_charlen_t l1);
    extern void FC_GLOBAL(sgges,SGGES)(char *JOBVSL, char*JOBVSRL, char * SORT, void * SELECTG, Int *N, float *A, Int *LDA, float *C, Int *ldc, Int *SDIM, float *ALPHAR, float *ALPHAI, float *betar,
            float *VSL, Int *LDVSL, float *VSR, Int *LDVR,  float *WORK, Int * LDWORK, Int *BWORK, Int *INFO);
    extern void FC_GLOBAL(sgehrd, SGEHRD)(Int *N, Int *ILO, Int *IHI, float *A, Int *LDA, float *TAU, float *WORK, Int *LDWORK, Int *Info);
    extern void FC_GLOBAL(sggbal, SGGBAL)(char* JOB, Int *N, float *A, Int *LDA, float *B, Int *LDB, Int *ILO, Int *IHI, float *LSCALE, float *RSCALE, float *WORK, Int *INFO, fortran_charlen_t l1);
    extern void FC_GLOBAL(sgghrd, SGGHRD)(char *COMPQ, char* COMPZ, Int *N, Int *ILO, Int *IHI, float *A, Int *LDA,
                    float *B, Int *LDB, float *Q, Int *LDQ, float *Z, Int *LDZ, Int *INFO, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(sgeqrf, SGEQRF)(Int *M, Int *N, float *A, Int *LDA, float *TAU, float *WORK, Int *LWORK, Int *INFO);
    extern void FC_GLOBAL(sormqr, SORMQR)(char *SIDE, char *TRANS, Int *M, Int *N, Int *K, float *A, Int *LDA, float *TAU, float *C, Int *LDC, float *WORK, Int *LWORK, Int *INFO, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(sorgqr, SORGQR)(Int *M, Int *N, Int *K, float *A, Int *LDA, float *TAU, float *WORK, Int *LWORK, Int *INFO);

    /*-----------------------------------------------------------------------------
     *  Error Handler
     *-----------------------------------------------------------------------------*/
    void FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C) ( void (*eh)(const char*,int));
    void FC_GLOBAL_(xerror_set_handler_f,XERROR_SET_HANDLER_F) ( void (*eh)(const char*,Int*, fortran_charlen_t));

    /*-----------------------------------------------------------------------------
     *  MEPACK Double
     *-----------------------------------------------------------------------------*/


    MEPACK_EXPORT void FC_GLOBAL_(dla_sort_gev,DLA_SORT_GEV)( Int * M,  double * A, Int * lDA, double *C, Int *ldc, double * Q, Int *ldq, double *Z, Int *ldz, Int *NB, double * work, Int *ldwork, Int *info );
    MEPACK_EXPORT void FC_GLOBAL_(dla_sort_ev,DLA_SORT_EV)( Int * M,  double * A, Int * lDA, double * Q, Int *ldq, Int *NB, double * work, Int *ldwork, Int *info );
    MEPACK_EXPORT void FC_GLOBAL_(sla_sort_ev,SLA_SORT_EV)( Int * M,  float * A, Int * lDA, float * Q, Int *ldq, Int *NB, float * work, Int *ldwork, Int *info );
    MEPACK_EXPORT void FC_GLOBAL_(sla_sort_gev,SLA_SORT_GEV)( Int * M,  float * A, Int * lDA, float *C, Int *ldc, float * Q, Int *ldq, float *Z, Int *ldz, Int *NB, float * work, Int *ldwork, Int *info );

    /*-----------------------------------------------------------------------------
     *  RECSY
     *-----------------------------------------------------------------------------*/
    extern void FC_GLOBAL(recsyct,RECSYCT)( Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double * E, Int *LDE, Int * INFO, double * MACHINE);
    extern void FC_GLOBAL(recsydt,RECSYDT)( Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double * E, Int *LDE, Int * INFO, double * MACHINE, double *work, Int *worksize);
    extern void FC_GLOBAL(recgsyl,RECGSYL)( Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double * E, Int *LDE, Int * INFO, double * MACHINE, double *work, Int *worksize);
    extern void FC_GLOBAL(recglyct,RECGLYCT)(Int *UPLO, double *SCALE, Int *N, double *A, Int *LDA, double *E, Int *LDE, double *X, Int *LDX, Int *info,  double *MACHINE, double *WORK, Int *LDWORK);
    extern void FC_GLOBAL(recglydt,RECGLYDT)(Int *UPLO, double *SCALE, Int *N, double *A, Int *LDA, double *E, Int *LDE, double *X, Int *LDX, Int *info,  double *MACHINE, double *WORK, Int *LDWORK);
    extern void FC_GLOBAL(recgcsy,RECGCSY)( Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double * E, Int *LDE, double *F, Int *LDF,  Int * INFO, double * MACHINE);
    extern void FC_GLOBAL(reclydt,RECLYDT)( Int * UPLOSIGN, double * SCALE, Int *M, double *A, Int *LDA, double * E, Int *LDE, Int * INFO, double * MACHINE, double *work, Int *worksize);
    extern void FC_GLOBAL(reclyct,RECLYCT)( Int * UPLOSIGN, double * SCALE, Int *M, double *A, Int *LDA, double * E, Int *LDE, Int * INFO, double * MACHINE);

    extern void FC_GLOBAL_(recsyct_p,RECSYCT_P)( Int *procs, Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double * E, Int *LDE, Int * INFO, double * MACHINE);
    extern void FC_GLOBAL_(recsydt_p,RECSYDT_P)( Int *proc, Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double * E, Int *LDE, Int * INFO, double * MACHINE, double *work, Int *worksize);
    extern void FC_GLOBAL_(recgsyl_p,RECGSYL_P)( Int *procs,  Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double * E, Int *LDE, Int * INFO, double * MACHINE, double *work, Int *worksize);
    extern void FC_GLOBAL_(recglyct_p,RECGLYCT_P)(Int *proc, Int *UPLO, double *SCALE, Int *N, double *A, Int *LDA, double *E, Int *LDE, double *X, Int *LDX, Int *info,  double *MACHINE, double *WORK, Int *LDWORK);
    extern void FC_GLOBAL_(recglydt_p,RECGLYDT_P)(Int *proc, Int *UPLO, double *SCALE, Int *N, double *A, Int *LDA, double *E, Int *LDE, double *X, Int *LDX, Int *info,  double *MACHINE, double *WORK, Int *LDWORK);
    extern void FC_GLOBAL_(recgcsy_p,RECGCSY_P)( Int *proc, Int * UPLOSIGN, double * SCALE, Int *M, Int *N, double *A, Int *LDA,double * B, Int * LDB, double *C, Int *LDC, double *D, Int *LDD, double * E, Int *LDE, double *F, Int *LDF,  Int * INFO, double * MACHINE);
    extern void FC_GLOBAL_(reclydt_p,RECLYDT_P)( Int *proc, Int * UPLOSIGN, double * SCALE, Int *M, double *A, Int *LDA, double * E, Int *LDE, Int * INFO, double * MACHINE, double *work, Int *worksize);
    extern void FC_GLOBAL_(reclyct_p,RECLYCT_P)( Int *proc, Int * UPLOSIGN, double * SCALE, Int *M, double *A, Int *LDA, double * E, Int *LDE, Int * INFO, double * MACHINE);

    extern void FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(double * MACHINE);

    /*-----------------------------------------------------------------------------
     *  SLICOT
     *-----------------------------------------------------------------------------*/
    extern void FC_GLOBAL(sg03ay,SG03AY)(char *TRANS, Int *N, double *A, Int *LDA, double *B, Int *LDB, double *X, Int *LDX, double *SCALE, Int *info, fortran_charlen_t l1);
    extern void FC_GLOBAL(sg03ax,SG03AX)(char *TRANS, Int *N, double *A, Int *LDA, double *B, Int *LDB, double *X, Int *LDX, double *SCALE, Int *info, fortran_charlen_t l1);
    extern void FC_GLOBAL(sg03ad,SG03AD)(char *DICO, char *JOB, char *FACT, char *TRANS, char *UPLO,
            Int *N, double *A, Int *LDA, double *B, Int *LDB, double *Q, Int *LDQ, double *Z, Int *LDZ,  double *X, Int *LDX,
            double *SCALE, double *SEP, double *FERR, double*ALPHAR, double *ALPHAI, double *BETA, Int *IWORK,
            double *DWORK, Int *LDWORK, Int *info, fortran_charlen_t l1, fortran_charlen_t l2, fortran_charlen_t l3, fortran_charlen_t l4, fortran_charlen_t l5);
    extern void FC_GLOBAL(sb03md,SB03MD)(char *DICO, char *JOB, char *FACT, char *TRANS,
            Int *N, double *A, Int *LDA, double *Q, Int *LDQ, double *X, Int *LDX,
            double *SCALE, double *SEP, double *FERR, double*ALPHAR, double *ALPHAI, Int *IWORK,
            double *DWORK, Int *LDWORK, Int *info, fortran_charlen_t l1, fortran_charlen_t l2, fortran_charlen_t l3, fortran_charlen_t l4);
    extern void FC_GLOBAL(sb04md,SB04MD)( Int *N, Int *M, double *A, Int *LDA, double *B, Int *LDB, double *C, Int *LDC, double *Z, Int *LDZ, Int *IWORK, double *DWORK, Int *LDWORK, Int *INFO);
    extern void FC_GLOBAL(sb04nd,SB04ND)( char *ABSCHU, char *ULA, char *ULB, Int *N, Int *M, double *A, Int *LDA, double *B, Int *LDB, double *C, Int *LDC, double *TOL, Int *IWORK, double *DWORK, Int *LDWORK, Int *INFO, fortran_charlen_t l1, fortran_charlen_t l2, fortran_charlen_t l3);
    extern void FC_GLOBAL(sb04qd,SB04QD)( Int *N, Int *M, double *A, Int *LDA, double *B, Int *LDB, double *C, Int *LDC, double *Z, Int *LDZ, Int *IWORK, double *DWORK, Int *LDWORK, Int *INFO);
    extern void FC_GLOBAL(sb04rd,SB04RD)( char *ABSCHU, char *ULA, char *ULB, Int *N, Int *M, double *A, Int *LDA, double *B, Int *LDB, double *C, Int *LDC, double *TOL, Int *IWORK, double *DWORK, Int *LDWORK, Int *INFO, fortran_charlen_t l1 , fortran_charlen_t l2, fortran_charlen_t l3);
    extern void FC_GLOBAL(sb04py,SB04PY)(char * TRANSA, char *TRANSB, Int *ISGN, Int *M, Int *N, double *A, Int *LDA, double *B, Int *LDB, double *X, Int *LDX, double *SCALE, double *dwork, Int *Info, fortran_charlen_t l1, fortran_charlen_t l2);
    extern void FC_GLOBAL(sb03my,SB03MY)( char *TRANA, Int * N, double *A, Int *LDA, double *C, Int *LDC, double *SCAL, Int *INFO, fortran_charlen_t l1);
    extern void FC_GLOBAL(sb03mx,SB03MX)( char *TRANA, Int * N, double *A, Int *LDA, double *C, Int *LDC, double *SCAL, double *WORK, Int *INFO, fortran_charlen_t l1);

    /* Algorithm 705  */
    extern void FC_GLOBAL(bkcon,BKCON)(Int *N, double *A, Int *LDA, double *B, Int *LDB, double *X, Int *LDX, Int *JOB, Int *info);
    extern void FC_GLOBAL(bkdis,BKDIS)(Int *N, double *A, Int *LDA, double *B, Int *LDB, double *X, Int *LDX, Int *JOB, Int *info);
    extern void FC_GLOBAL(bkhs2,BKHS2)(Int *M, Int *N, double *A, Int *LDA, double *B, Int *LDB,  double *C, Int *LDC, double *D, Int *LDD, double *X, Int *LDX, double *WORK, Int *iwork, Int *JOB, Int *info);

    /*-----------------------------------------------------------------------------
     *  LAPACK
     *-----------------------------------------------------------------------------*/
    extern void FC_GLOBAL(dtrsyl,DTRSYL)(char * TRANSA, char *TRANSB, Int *ISGN, Int *M, Int *N, double *A, Int *LDA, double *B, Int *LDB, double *X, Int *LDX, double *SCALE, Int *Info);
    extern void FC_GLOBAL(strsyl,STRSYL)(char * TRANSA, char *TRANSB, Int *ISGN, Int *M, Int *N, float *A, Int *LDA, float *B, Int *LDB, float *X, Int *LDX, float *SCALE, Int *Info);




    typedef void (*context_init_t)(void *_ctx, Int M, Int N, Int ISolver, Int MB, Int NB, Int BIGMB, const char *TA, const char *TB, double sgn1, double sgn2);
    typedef void (*context_setup_rhs_t)(void *_ctx);
    typedef void (*context_setup_problem_t)(void *_ctx, Int iseed[4], Int align, Int changerole);
    typedef void (*context_iter_prepare_t)(void * _ctx);
    typedef void (*context_iter_solve_t)(void *_ctx);
    typedef void (*context_check_t)(void *_ctx, double *forward, double *rel);
    typedef void (*context_clean_t)(void *_ctx);

    typedef struct _benchmark_context_t {
        void *ctx;
        context_init_t init;
        context_setup_rhs_t setup_rhs;
        context_setup_problem_t setup_problem;
        context_iter_prepare_t iter_prepare;
        context_iter_solve_t iter_solve;
        context_check_t check;
        context_clean_t clean;
    } benchmark_context_t;

    void benchmark_context_free(benchmark_context_t *bctx);
    void benchmark_context_init(benchmark_context_t *ctx, Int M, Int N, Int ISolver, Int MB, Int NB, Int BIGMB, const char *TA, const char *TB, double sgn1, double sgn2);
    void benchmark_context_setup_rhs(benchmark_context_t  *ctx);
    void benchmark_context_setup_problem(benchmark_context_t *ctx, Int iseed[4], Int align, Int changerole);
    void benchmark_context_iter_prepare(benchmark_context_t *ctx);
    void benchmark_context_iter_solve(benchmark_context_t *ctx);
    void benchmark_context_check(benchmark_context_t *ctx, double *forward, double *rel);
    void benchmark_context_clean(benchmark_context_t *ctx);

#ifdef __cplusplus
};
#endif
