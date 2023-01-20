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

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "benchmark.h"

double benchmark_check_X_double(Int M, Int N, double *X1, Int LDX1, double *X2, Int LDX2)
{
    Int i, j;
    double e1 = 0, e2 = 0;
    double dummy = 0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            X1[i*LDX1+j] -= X2[i*LDX2+j];
        }
    }
    e1 = FC_GLOBAL(dlange,DLANGE)("1", &M, &N, X1, &LDX1, &dummy, 1) ;
    e2 = FC_GLOBAL(dlange,DLANGE)("1", &M, &N, X2, &LDX2, &dummy, 1) ;
    return e1/e2;
}

double FC_GLOBAL(benchmark_check_x_double_f,BENCHMARK_CHECK_X_DOUBLE_F)
    (Int *M, Int *N, double *X1, Int *LDX1, double *X2, Int *LDX2)
{
    return benchmark_check_X_double(*M, *N, X1, *LDX1, X2, *LDX2);
}

float benchmark_check_X_float(Int M, Int N, float *X1, Int LDX1, float *X2, Int LDX2)
{
    Int i, j;
    float nrm1, nrm2 ;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            X1[i*LDX1+j] -= X2[i*LDX2+j];
            /* printf("%d %d  = %g\n", j, i, X1[i*LDX1+j]); */
        }
    }
    nrm1 = FC_GLOBAL(slange,SLANGE)("1", &M, &N, X1, &M, NULL, 1);
    nrm2 = FC_GLOBAL(slange,SLANGE)("1", &M, &N, X2, &M, NULL, 1);

    return nrm1/nrm2;
}


void benchmark_rhs_sylv_double(char *TRANSA, char *TRANSB, double sgn,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY)
{
        double alpha, beta;

        alpha=1;
        beta =0;
        FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, X, &LDX, &beta, Y, &LDY, 1, 1);
        beta = 1;
        alpha = sgn;
        FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &N, &N, &alpha , X, &LDX, B, &LDB, &beta, Y, &LDY, 1, 1);
        return;

}

void benchmark_rhs_sylv_float(char *TRANSA, char *TRANSB, float sgn,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY)
{
        float alpha, beta;

        alpha=1;
        beta =0;
        FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, X, &LDX, &beta, Y, &LDY, 1, 1);
        beta = 1;
        alpha = sgn;
        FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &N, &N, &alpha , X, &LDX, B, &LDB, &beta, Y, &LDY, 1, 1);
        return;

}


void benchmark_rhs_sylv2_double(char *TRANSA, char *TRANSB, double sgn,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY)
{
        double alpha, beta;
        double *Work;

        Work = malloc(M*N*sizeof(double));
        alpha=1;
        beta =0;
        FC_GLOBAL(dlacpy,DLACPY)("All", &M, &N, X, &M, Y, &LDY, 1);
        FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
        alpha = 1;
        beta = sgn;
        FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &N, &N, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
        free(Work);
        return;

}

void benchmark_rhs_sylv2_float(char *TRANSA, char *TRANSB, float sgn,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY)
{
        float alpha, beta;
        float *Work;

        Work = malloc(M*N*sizeof(float));
        alpha=1;
        beta =0;
        FC_GLOBAL(slacpy,SLACPY)("All", &M, &N, X, &M, Y, &LDY, 1);
        FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
        alpha = 1;
        beta = sgn;
        FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &N, &N, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
        free(Work);
        return;

}



void benchmark_rhs_gsylv_double(char *TRANSA, char *TRANSB, double sgn,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double * C, Int LDC, double *D, Int LDD,
        double *X, Int LDX, double *Y, Int LDY)
{
        double alpha, beta;
        double *Work;

        Work = malloc(M*N*sizeof(double));
        alpha=1;
        beta =0;
        FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
        FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &N, &N, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);

        FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &N, &M, &alpha , C, &LDC, X, &LDX, &beta, Work, &M, 1, 1);
        beta = 1;
        alpha = sgn;
        FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &N, &N, &alpha , Work, &M, D, &LDD, &beta, Y, &LDY, 1, 1);
        free(Work);
        return;

}

void FC_GLOBAL(benchmark_rhs_gsylv_double_f,BENCHMARK_RHS_GSYLV_DOUBLE_F)
        (char *TRANSA, char *TRANSB, double *sgn,  Int *M, Int *N, double * A, Int *LDA, double * B, Int *LDB, double * C, Int *LDC, double *D, Int *LDD,
         double *X, Int *LDX, double *Y, Int *LDY, fortran_charlen_t l1, fortran_charlen_t l2 )
{
    (void ) l1;
    (void ) l2;
    char ta[2] = { TRANSA[0], 0};
    char tb[2] = { TRANSB[0], 0};

    benchmark_rhs_gsylv_double(ta, tb, *sgn,  *M, *N, A, *LDA, B, *LDB, C, *LDC, D, *LDD, X, *LDX, Y, *LDY);
}



void benchmark_rhs_gsylv_float(char *TRANSA, char *TRANSB, float sgn,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float * C, Int LDC, float *D, Int LDD,
        float *X, Int LDX, float *Y, Int LDY)
{
        float alpha, beta;
        float *Work;

        Work = malloc(M*N*sizeof(float));
        alpha=1;
        beta =0;
        FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
        FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &N, &N, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
        FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &N, &M, &alpha , C, &LDC, X, &LDX, &beta, Work, &M, 1, 1);
        beta = 1;
        alpha = sgn;
        FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &N, &N, &alpha , Work, &M, D, &LDD, &beta, Y, &LDY, 1, 1);
        free(Work);
        return;

}


void benchmark_rhs_ggcsylv_double(char *TRANSA, char *TRANSB, double sgn1, double sgn2,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double * C, Int LDC, double *D, Int LDD,
        double *R, Int LDR, double *L, Int LDL, double *Y1, Int LDY1, double *Y2, Int LDY2)
{
        double alpha, beta;
        alpha=1;
        beta =0;
        FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, R, &LDR, &beta, Y1, &LDY1, 1, 1);
        beta = 1;
        alpha = sgn1;
        FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &N, &N, &alpha,  L, &LDL, B, &LDB, &beta, Y1, &LDY1, 1, 1);

        alpha = 1;
        beta =0;
        FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &N, &M, &alpha , C, &LDC, R, &LDR, &beta, Y2, &LDY2, 1, 1);
        beta = 1;
        alpha = sgn2;
        FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &N, &N, &alpha,  L, &LDL, D, &LDD, &beta, Y2, &LDY2, 1, 1);

        return;

}

void benchmark_rhs_ggcsylv_dual_double(char *TRANSA, char *TRANSB, double sgn1, double sgn2,  Int M, Int N, double * A, Int LDA, double * B, Int LDB, double * C, Int LDC, double *D, Int LDD,
        double *R, Int LDR, double *L, Int LDL, double *Y1, Int LDY1, double *Y2, Int LDY2)
{
        double alpha, beta;
        char TA[2] = {0 , 0};
        char TB[2] = {0 , 0};

        if ( tolower (TRANSA[0]) == 't' ) {
            TA [0] = 'N';
        } else {
            TA [0] = 'T';
        }

        if ( tolower (TRANSB[0]) == 't' ) {
            TB [0] = 'N';
        } else {
            TB [0] = 'T';
        }

        alpha = 1;
        beta = 0;
        FC_GLOBAL(dgemm,DGEMM)(TA, "N", &M, &N, &M, &alpha , A, &LDA, R, &LDR, &beta, Y1, &LDY1, 1, 1);

        beta = 1;
        alpha = 1;
        FC_GLOBAL(dgemm,DGEMM)(TA, "N", &M, &N, &M, &alpha,  C, &LDC, L, &LDL, &beta, Y1, &LDY1, 1, 1);

        beta = 0;
        alpha = sgn1;
        FC_GLOBAL(dgemm,DGEMM)("N", TB, &M, &N, &N, &alpha,  R, &LDR, B, &LDB, &beta, Y2, &LDY2, 1, 1);

        beta = 1;
        alpha = sgn2;
        FC_GLOBAL(dgemm,DGEMM)("N", TB, &M, &N, &N, &alpha,  L, &LDL, D, &LDD, &beta, Y2, &LDY2, 1, 1);

        return;

}

void benchmark_rhs_ggcsylv_float(char *TRANSA, char *TRANSB, float sgn1, float sgn2,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float * C, Int LDC, float *D, Int LDD,
        float *R, Int LDR, float *L, Int LDL, float *Y1, Int LDY1, float *Y2, Int LDY2)
{
        float alpha, beta;
        alpha=1;
        beta =0;
        FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &N, &M, &alpha , A, &LDA, R, &LDR, &beta, Y1, &LDY1, 1, 1);
        beta = 1;
        alpha = sgn1;
        FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &N, &N, &alpha,  L, &LDL, B, &LDB, &beta, Y1, &LDY1, 1, 1);

        alpha = 1;
        beta =0;
        FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &N, &M, &alpha , C, &LDC, R, &LDR, &beta, Y2, &LDY2, 1, 1);
        beta = 1;
        alpha = sgn2;
        FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &N, &N, &alpha,  L, &LDL, D, &LDD, &beta, Y2, &LDY2, 1, 1);

        return;

}


void benchmark_rhs_ggcsylv_dual_float(char *TRANSA, char *TRANSB, float sgn1, float sgn2,  Int M, Int N, float * A, Int LDA, float * B, Int LDB, float * C, Int LDC, float *D, Int LDD,
        float *R, Int LDR, float *L, Int LDL, float *Y1, Int LDY1, float *Y2, Int LDY2)
{
        float alpha, beta;
        char TA[2] = {0 , 0};
        char TB[2] = {0 , 0};

        if ( tolower (TRANSA[0]) == 't' ) {
            TA [0] = 'N';
        } else {
            TA [0] = 'T';
        }

        if ( tolower (TRANSB[0]) == 't' ) {
            TB [0] = 'N';
        } else {
            TB [0] = 'T';
        }

        alpha = 1;
        beta = 0;
        FC_GLOBAL(sgemm,SGEMM)(TA, "N", &M, &N, &M, &alpha , A, &LDA, R, &LDR, &beta, Y1, &LDY1, 1, 1);

        beta = 1;
        alpha = 1;
        FC_GLOBAL(sgemm,SGEMM)(TA, "N", &M, &N, &M, &alpha,  C, &LDC, L, &LDL, &beta, Y1, &LDY1, 1, 1);

        beta = 0;
        alpha = sgn1;
        FC_GLOBAL(sgemm,SGEMM)("N", TB, &M, &N, &N, &alpha,  R, &LDR, B, &LDB, &beta, Y2, &LDY2, 1, 1);

        beta = 1;
        alpha = sgn2;
        FC_GLOBAL(sgemm,SGEMM)("N", TB, &M, &N, &N, &alpha,  L, &LDL, D, &LDD, &beta, Y2, &LDY2, 1, 1);

        return;

}



void benchmark_rhs_glyap_double(char *TRANSA, Int M, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY)
{
    double alpha, beta;
    double *Work;
    char TRANSB[2];

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';

    Work = malloc(M*M*sizeof(double));
    alpha = 1;
    beta  = 0;
    FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
    FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
    FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &M, &M, &alpha , B, &LDB, X, &LDX, &beta, Work, &M, 1, 1);
    beta = 1;
    alpha = 1;
    FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, A, &LDA, &beta, Y, &LDY, 1, 1);
    free(Work);
    return;

}

void benchmark_rhs_glyap_float(char *TRANSA, Int M, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY)
{
    float alpha, beta;
    float *Work;
    char TRANSB[2];

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';

    Work = malloc(M*M*sizeof(float));
    alpha = 1;
    beta  = 0;
    FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
    FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
    FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &M, &M, &alpha , B, &LDB, X, &LDX, &beta, Work, &M, 1, 1);
    beta = 1;
    alpha = 1;
    FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, A, &LDA, &beta, Y, &LDY, 1, 1);
    free(Work);
    return;

}


void benchmark_rhs_lyap_double(char *TRANSA, Int M, double * A, Int LDA, double *X, Int LDX, double *Y, Int LDY)
{
    double alpha, beta;
    char TRANSB[2];

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';

    alpha = 1;
    beta  = 0;
    FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, Y , &LDY, 1, 1);
    beta = 1;
    FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &M, &M, &alpha , X, &LDX, A, &LDA, &beta, Y, &LDY, 1, 1);
    return;

}

void benchmark_rhs_lyap_float(char *TRANSA, Int M, float * A, Int LDA, float *X, Int LDX, float *Y, Int LDY)
{
    float alpha, beta;
    char TRANSB[2];

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';

    alpha = 1;
    beta  = 0;
    FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, Y , &LDY, 1, 1);
    beta = 1;
    FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &M, &M, &alpha , X, &LDX, A, &LDA, &beta, Y, &LDY, 1, 1);
    return;

}

void benchmark_rhs_stein_double(char *TRANSA, Int M, double * A, Int LDA, double *X, Int LDX, double *Y, Int LDY)
{
    double alpha, beta;
    char TRANSB[2];
    double *work;
    Int i, j;

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';
    work = malloc(sizeof(double) *M*M);

    alpha = 1;
    beta  = 0;
    FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, work, &M, 1, 1);
    FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &M, &M, &alpha , work, &M, A, &LDA, &beta, Y, &LDY, 1, 1);
    for (j = 0; j < M; j++) {
        for (i = 0; i < M; i++) {
            Y[i+j*LDY] -= X[i+j*LDX];
        }
    }

    free(work);
    return;


}

void benchmark_rhs_stein_float(char *TRANSA, Int M, float * A, Int LDA, float *X, Int LDX, float *Y, Int LDY)
{
    float alpha, beta;
    char TRANSB[2];
    float *work;
    Int i, j;

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';
    work = malloc(sizeof(float) *M*M);

    alpha = 1;
    beta  = 0;
    FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, work, &M, 1, 1);
    FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &M, &M, &alpha , work, &M, A, &LDA, &beta, Y, &LDY, 1, 1);
    for (j = 0; j < M; j++) {
        for (i = 0; i < M; i++) {
            Y[i+j*LDY] -= X[i+j*LDX];
        }
    }

    free(work);
    return;


}



void benchmark_rhs_gstein_double(char *TRANSA, Int M, double * A, Int LDA, double * B, Int LDB, double *X, Int LDX, double *Y, Int LDY)
{
    double alpha, beta;
    double *Work;
    char TRANSB[2];

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';

    Work = malloc(M*M*sizeof(double));
    alpha = 1;
    beta  = 0;
    FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
    FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, A, &LDA, &beta, Y, &LDY, 1, 1);
    FC_GLOBAL(dgemm,DGEMM)(TRANSA, "N", &M, &M, &M, &alpha , B, &LDB, X, &LDX, &beta, Work, &M, 1, 1);
    beta = 1;
    alpha = -1;
    FC_GLOBAL(dgemm,DGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
    free(Work);
    return;

}


void benchmark_rhs_gstein_float(char *TRANSA, Int M, float * A, Int LDA, float * B, Int LDB, float *X, Int LDX, float *Y, Int LDY)
{
    float alpha, beta;
    float *Work;
    char TRANSB[2];

    if ( TRANSA[0] == 'N')
        TRANSB[0] = 'T';
    else
        TRANSB[0] = 'N';

    Work = malloc(M*M*sizeof(float));
    alpha = 1;
    beta  = 0;
    FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &M, &M, &alpha , A, &LDA, X, &LDX, &beta, Work, &M, 1, 1);
    FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, A, &LDA, &beta, Y, &LDY, 1, 1);
    FC_GLOBAL(sgemm,SGEMM)(TRANSA, "N", &M, &M, &M, &alpha , B, &LDB, X, &LDX, &beta, Work, &M, 1, 1);
    beta = 1;
    alpha = -1;
    FC_GLOBAL(sgemm,SGEMM)("N", TRANSB, &M, &M, &M, &alpha , Work, &M, B, &LDB, &beta, Y, &LDY, 1, 1);
    free(Work);
    return;

}


