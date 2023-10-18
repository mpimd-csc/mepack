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
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include <time.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include "benchmark.h"
#include "mepack.h"
#include "mepack_internal.h"

#include "cscutils/table.h"

#ifdef SINGLE_PRECISION
#define STRINGIFY(x) #x
#define MEPACK_PRECISION_PREFIX(X) STRINGIFY( SLA_ ## X)
#define MEPACK_PREFIX(x) mepack_single_ ## x
#define FLOAT float
#else
#define MEPACK_PREFIX(x) mepack_double_ ## x
#define FLOAT double
#define STRINGIFY(x) #x
#define MEPACK_PRECISION_PREFIX(X) STRINGIFY( DLA_ ## X)
#endif


typedef struct _context_sylv_t {
    Int M, N;
    Int MB, NB, BIGMB;
    FLOAT sgn;
    FLOAT scale;
    FLOAT *A;  Int LDA;
    FLOAT *B;  Int LDB;
    FLOAT *X;  Int LDX;
    FLOAT *Xorig; Int LDXorig;
    FLOAT *RHS_X; Int LDRHS_X;
    ssize_t mem;
    FLOAT *work;
    char * TRANSA;
    char * TRANSB;
    int isolver;
#ifdef RECSY
    double MACHINE_RECSY[12];
    Int type;
#endif
    int changerole;
} context_sylv_t;


int is_sym(void) {
    return 0;
}


static void context_sylv_init(void *_ctx, Int M, Int N, Int ISolver, Int MB, Int NB, Int BIGMB, const char *TA, const char *TB, double sgn1, double sgn2){
    context_sylv_t *ctx = _ctx;
    Int i;
    (void) sgn2;

    MEPACK_PREFIX(trsylv_blocksize_2stage_set)(BIGMB);
    MEPACK_PREFIX(trsylv_blocksize_mb_set)(MB);
    MEPACK_PREFIX(trsylv_blocksize_nb_set)(NB);

    if ( ISolver < 19 ) {
        if ( ISolver < 7 )  {
            mepack_trsylv_isolver_set(ISolver);
        } else if ( ISolver > 6 && ISolver < 13 )  {
            mepack_trsylv_isolver_set(ISolver-6);
        } else if ( ISolver > 12 && ISolver < 19 )  {
            mepack_trsylv_isolver_set(ISolver-12);
        } else {
            mepack_trsylv_isolver_set(1);
        }
    } else {
        mepack_trsylv_isolver_set(1);
    }

    ctx->A = (FLOAT *) malloc(sizeof(FLOAT) * (M*M));
    ctx->B = (FLOAT *) malloc(sizeof(FLOAT) * (N*N));
    ctx->X = (FLOAT *) malloc(sizeof(FLOAT) * (M*N));


#pragma omp parallel for schedule(static,1)
    for (i = 0; i < M*M; i+=512) { ctx->A[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
    for (i = 0; i < N*N; i+=512) { ctx->B[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
    for (i = 0; i < M*N; i+=512) { ctx->X[i] = 0.0; }

    ctx->Xorig = (FLOAT *) malloc(sizeof(FLOAT) * (M*N));
    ctx->RHS_X = (FLOAT *) malloc(sizeof(FLOAT) * (M*N));

    if ( ISolver < 7 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L3), M, N);
    } else if ( ISolver > 6 && ISolver < 13 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_DAG), M, N);
    } else if ( ISolver > 12 && ISolver < 19 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L3_2S), M, N);
    } else if ( ISolver == 19 ) {
        if ( M <=32 && N <= 32)
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_LOCAL_COPY_32), M, N);
        else if ( M <=64 && N <= 64)
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_LOCAL_COPY_64), M, N);
        else if ( M <=96 && N <= 96)
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_LOCAL_COPY_96), M, N);
        else if ( M <=128 && N <= 128)
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_LOCAL_COPY_128), M, N);
        else
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_REORDER), M, N);
    } else if ( ISolver == 20 ) {
        if ( M <=128 && N <= 128)
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_LOCAL_COPY), M, N);
        else
            ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_REORDER), M, N);
    } else if ( ISolver == 21)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_REORDER), M, N);
    else if ( ISolver == 22)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2), M, N);
    else if ( ISolver == 23)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_UNOPT), M, N);
    else if ( ISolver == 24)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_RECURSIVE), M, N);
    else
        ctx->mem = 2*M*N;

    if ( ISolver == 19 || ISolver == 20 ) {
        ctx->mem = MAX(mepack_memory(MEPACK_PRECISION_PREFIX(TRSYLV_L2_REORDER), M, N), ctx->mem);
    }
    ctx->work = (FLOAT *) malloc(sizeof(FLOAT) * (ctx->mem));

    ctx->M = M;
    ctx->N = N;
    ctx->MB = MB;
    ctx->NB = NB;
    ctx->BIGMB = BIGMB;
    ctx->sgn = sgn1;
    ctx->scale = 1;
    ctx->LDA = M;
    ctx->LDB = N;
    ctx->LDX = M;
    ctx->LDXorig = M;
    ctx->LDRHS_X = M;
    ctx->isolver = ISolver;

    ctx->TRANSA=strdup(TA);
    ctx->TRANSB=strdup(TB);
#ifdef RECSY
    ctx->type = benchmark_sylv_setup_to_type(TA, TB, sgn1);
    FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(ctx->MACHINE_RECSY);

#endif

}

static void context_sylv_setup_rhs(void *_ctx){
    context_sylv_t * ctx = _ctx;

    FLOAT alpha = 1;
    FLOAT beta = 1;

#ifdef SINGLE_PRECISION
    FC_GLOBAL_(slaset,SLASET)("All", &ctx->M, &ctx->N, &alpha, &beta, ctx->Xorig, &ctx->LDXorig, 1);
#else
    FC_GLOBAL_(dlaset,DLASET)("All", &ctx->M, &ctx->N, &alpha, &beta, ctx->Xorig, &ctx->LDXorig, 1);
#endif
}

static void context_sylv_setup_problem(void *_ctx, Int iseed[4], Int align, Int changerole)
{
    context_sylv_t * ctx = _ctx;
    (void) changerole;
#ifdef SINGLE_PRECISION

    benchmark_random_evp_float(ctx->M, iseed, ctx->A, NULL, NULL, align);
    benchmark_random_evp_float(ctx->N, iseed, ctx->B, NULL, NULL, align);

    benchmark_rhs_sylv_float(ctx->TRANSA, ctx->TRANSB, ctx->sgn, ctx->M, ctx->N,
            ctx->A, ctx->LDA, ctx->B, ctx->LDB,
            ctx->Xorig, ctx->LDXorig, ctx->RHS_X, ctx->LDRHS_X);
#else
    benchmark_random_evp_double(ctx->M, iseed, ctx->A, NULL, NULL, align);
    benchmark_random_evp_double(ctx->N, iseed, ctx->B, NULL, NULL, align);

    benchmark_rhs_sylv_double(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
            ctx->A, ctx->LDA, ctx->B, ctx->LDB,
            ctx->Xorig, ctx->LDXorig, ctx->RHS_X, ctx->LDRHS_X);

#endif
}

static void context_sylv_iter_prepare(void * _ctx) {
    context_sylv_t * ctx = _ctx;
#ifdef SINGLE_PRECISION
    FC_GLOBAL_(slacpy,SLACPY)("All", &ctx->M, &ctx->N, ctx->RHS_X, &ctx->LDRHS_X, ctx->X, &ctx->LDX, 1);
#else
    FC_GLOBAL_(dlacpy,DLACPY)("All", &ctx->M, &ctx->N, ctx->RHS_X, &ctx->LDRHS_X, ctx->X, &ctx->LDX, 1);
#endif
}

static void context_sylv_iter_solve(void *_ctx){
    context_sylv_t * ctx = _ctx;
    int info= 0;

    if ( ctx->isolver < 7 ) {
        MEPACK_PREFIX(trsylv_level3)(ctx->TRANSA, ctx->TRANSB, ctx->sgn, ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver < 12 ) {
        MEPACK_PREFIX(trsylv_dag)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver < 19 ) {
        MEPACK_PREFIX(trsylv_level3_2stage)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 19 ){
        if ( ctx -> M <= 32 && ctx->N <= 32 )
            MEPACK_PREFIX(trsylv_level2_local_copy_32)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
        else if ( ctx->M <= 64 && ctx->N <= 64 )
            MEPACK_PREFIX(trsylv_level2_local_copy_64)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
        else if ( ctx->M <= 96 && ctx->N <= 96 )
            MEPACK_PREFIX(trsylv_level2_local_copy_96)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
        else if ( ctx->M <= 128 && ctx->N <= 128 )
            MEPACK_PREFIX(trsylv_level2_local_copy_128)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
        else
            MEPACK_PREFIX(trsylv_level2_reorder)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 20 ) {
        if ( ctx->M <= 128 && ctx->N <= 128) {
            MEPACK_PREFIX(trsylv_level2_local_copy)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
        } else {
            MEPACK_PREFIX(trsylv_level2_reorder)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
        }
    } else if ( ctx->isolver == 21 ) {
        MEPACK_PREFIX(trsylv_level2_reorder)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 22 ) {
        MEPACK_PREFIX(trsylv_level2)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 23 ) {
        MEPACK_PREFIX(trsylv_level2_unopt)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 24){
        MEPACK_PREFIX(trsylv_recursive)(ctx->TRANSA, ctx->TRANSB, ctx->sgn,  ctx->M, ctx->N,
                ctx->A, ctx->LDA, ctx->B, ctx->LDB, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    }
#if defined(RECSY) && !defined (SINGLE_PRECISION)
    else if ( ctx->isolver == 25 ) {
        FC_GLOBAL_(recsyct,RECSYCT)(&ctx->type, &ctx->scale, &ctx->M, &ctx->N, ctx->A, &ctx->LDA, ctx->B, &ctx->LDB, ctx->X, &ctx->LDX,  &info, ctx->MACHINE_RECSY);
    } else if ( ctx->isolver == 26 ) {
        Int proc = omp_get_num_procs();
        FC_GLOBAL_(recsyct_p,RECSYCT_P)(&proc, &ctx->type, &ctx->scale, &ctx->M, &ctx->N, ctx->A, &ctx->LDA, ctx->B, &ctx->LDB, ctx->X, &ctx->LDX, &info, ctx->MACHINE_RECSY);
    }
#endif
    else if ( ctx->isolver == 29 ) {
        Int isgn;
        if ( ctx->sgn < 0 )
            isgn = -1;
        else
            isgn = 1;
    #ifdef SINGLE_PRECISION
        FC_GLOBAL_(strsyl,STRSYL)(ctx->TRANSA, ctx->TRANSB, &isgn, &ctx->M, &ctx->N, ctx->A, &ctx->LDA, ctx->B, &ctx->LDB,ctx->X,  &ctx->LDX, &ctx->scale, &info);
    #else
        FC_GLOBAL_(dtrsyl,DTRSYL)(ctx->TRANSA, ctx->TRANSB, &isgn, &ctx->M, &ctx->N, ctx->A, &ctx->LDA, ctx->B, &ctx->LDB,ctx->X,  &ctx->LDX, &ctx->scale, &info);
#endif

    }
}

static void context_sylv_check(void *_ctx, double *forward, double *rel)
{
    context_sylv_t * ctx = _ctx;
    *rel += MEPACK_PREFIX(residual_sylv)(ctx->TRANSA, ctx->TRANSB, ctx->sgn, ctx->M, ctx->N,
            ctx->A, ctx->LDA, ctx->B, ctx->LDB,
            ctx->X, ctx->LDX, ctx->RHS_X, ctx->LDRHS_X, ctx->scale);
#ifdef SINGLE_PRECISION
    *forward += benchmark_check_X_float(ctx->M,ctx->N,ctx->X, ctx->LDX, ctx->Xorig, ctx->LDXorig);
#else
    *forward += benchmark_check_X_double(ctx->M,ctx->N,ctx->X, ctx->LDX, ctx->Xorig, ctx->LDXorig);
#endif
}

static void context_sylv_clean(void *_ctx)
{
    context_sylv_t * ctx = _ctx;

    free(ctx->A);
    free(ctx->B);
    free(ctx->X);
    free(ctx->Xorig);
    free(ctx->RHS_X);
    free(ctx->work);
    free(ctx->TRANSA);
    free(ctx->TRANSB);
}

benchmark_context_t *benchmark_create(void)
{
    benchmark_context_t *bctx = malloc(sizeof(benchmark_context_t));
    context_sylv_t * ctx = malloc(sizeof(context_sylv_t));
    memset(ctx, 0, sizeof(context_sylv_t));
    bctx->ctx = ctx;
    bctx->init = context_sylv_init;
    bctx->setup_rhs = context_sylv_setup_rhs;
    bctx->setup_problem = context_sylv_setup_problem;
    bctx->iter_prepare = context_sylv_iter_prepare;
    bctx->iter_solve = context_sylv_iter_solve;
    bctx->check = context_sylv_check;
    bctx->clean = context_sylv_clean;
    return bctx;
}

const char *equation_name(void) {
    return "Triangular Generalized Sylvester Equation";
}

int have_solver(int is){
    if ( is > 0 && is <=29){
#ifndef RECSY
        if(is == 25 || is == 26 ) return 0;
#endif
        if (is == 27 || is == 28 ) return 0;
#ifdef SINGLE_PRECISION
        if(is == 25 || is == 26 ) return 0;
#endif
        return 1;
    }
    return 0;

}
