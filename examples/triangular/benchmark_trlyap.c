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


typedef struct _context_lyap_t {
    Int M ;
    Int MB, BIGMB;
    FLOAT sgn;
    FLOAT scale;
    FLOAT *A;  Int LDA;
    FLOAT *X;  Int LDX;
    FLOAT *Xorig; Int LDXorig;
    FLOAT *RHS_X; Int LDRHS_X;
    ssize_t mem;
    FLOAT *work;
    char * TRANSA;
    int isolver;
#ifdef RECSY
    double MACHINE_RECSY[12];
    Int type;
#endif
} context_lyap_t;


int is_sym(void) {
    return 1;
}


static void context_lyap_init(void *_ctx, Int M, Int N, Int ISolver, Int MB, Int NB, Int BIGMB, const char *TA, const char *TB, double sgn1, double sgn2){
    context_lyap_t *ctx = _ctx;
    Int i;
    (void) sgn1;
    (void) sgn2;
    (void) TB;
    (void) N;
    (void) NB;

    MEPACK_PREFIX(trlyap_blocksize_2stage_set)(BIGMB);
    MEPACK_PREFIX(trlyap_blocksize_set)(MB);

    if ( ISolver < 19 ) {
        if ( ISolver < 7 )  {
            mepack_trlyap_isolver_set(ISolver);
        } else if ( ISolver > 6 && ISolver < 13 )  {
            mepack_trlyap_isolver_set(ISolver-6);
        } else if ( ISolver > 12 && ISolver < 19 )  {
            mepack_trlyap_isolver_set(ISolver-12);
        } else {
            mepack_trlyap_isolver_set(1);
        }
    } else {
        mepack_trlyap_isolver_set(1);
    }

    ctx->A = (FLOAT *) malloc(sizeof(FLOAT) * (M*M));
    ctx->X = (FLOAT *) malloc(sizeof(FLOAT) * (M*M));


#pragma omp parallel for schedule(static,1)
    for (i = 0; i < M*M; i+=512) { ctx->A[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
    for (i = 0; i < M*M; i+=512) { ctx->X[i] = 0.0; }

    ctx->Xorig = (FLOAT *) malloc(sizeof(FLOAT) * (M*M));
    ctx->RHS_X = (FLOAT *) malloc(sizeof(FLOAT) * (M*M));

    if ( ISolver < 7 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L3), M, M);
    } else if ( ISolver > 6 && ISolver < 13 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_DAG), M, M);
    } else if ( ISolver > 12 && ISolver < 19 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L3_2S), M, M);
    } else if ( ISolver == 19 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L2_OPT), M, M);
    } else if ( ISolver == 20 ) {
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L2), M, M);
    } else if ( ISolver == 21)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L2), M, M);
    else if ( ISolver == 22)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L2), M, M);
    else if ( ISolver == 23)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_L2), M, M);
    else if ( ISolver == 24)
        ctx->mem = mepack_memory(MEPACK_PRECISION_PREFIX(TRLYAP_RECURSIVE), M, M);
    else
        ctx->mem = 3*M*M;

    ctx->work = (FLOAT *) malloc(sizeof(FLOAT) * (ctx->mem));

    ctx->M = M;
    ctx->MB = MB;
    ctx->BIGMB = BIGMB;
    ctx->scale = 1;
    ctx->LDA = M;
    ctx->LDX = M;
    ctx->LDXorig = M;
    ctx->LDRHS_X = M;
    ctx->isolver = ISolver;

    ctx->TRANSA=strdup(TA);
#ifdef RECSY
    ctx->type = benchmark_glyap_setup_to_type(TA);
    FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(ctx->MACHINE_RECSY);

#endif

}

static void context_lyap_setup_rhs(void *_ctx){
    context_lyap_t * ctx = _ctx;

    FLOAT alpha = 1;
    FLOAT beta = 1;

#ifdef SINGLE_PRECISION
    FC_GLOBAL_(slaset,SLASET)("All", &ctx->M, &ctx->M, &alpha, &beta, ctx->Xorig, &ctx->LDXorig, 1);
#else
    FC_GLOBAL_(dlaset,DLASET)("All", &ctx->M, &ctx->M, &alpha, &beta, ctx->Xorig, &ctx->LDXorig, 1);
#endif
}

static void context_lyap_setup_problem(void *_ctx, Int iseed[4], Int align, Int changerole)
{
    context_lyap_t * ctx = _ctx;
    (void) changerole;
#ifdef SINGLE_PRECISION

    benchmark_random_evp_float(ctx->M, iseed, ctx->A, NULL, NULL, align);

    benchmark_rhs_lyap_float(ctx->TRANSA, ctx->M,
            ctx->A, ctx->LDA,
            ctx->Xorig, ctx->LDXorig, ctx->RHS_X, ctx->LDRHS_X);
#else
    benchmark_random_evp_double(ctx->M, iseed, ctx->A, NULL, NULL, align);

    benchmark_rhs_lyap_double(ctx->TRANSA, ctx->M,
            ctx->A, ctx->LDA,
            ctx->Xorig, ctx->LDXorig, ctx->RHS_X, ctx->LDRHS_X);

#endif
}

static void context_lyap_iter_prepare(void * _ctx) {
    context_lyap_t * ctx = _ctx;
#ifdef SINGLE_PRECISION
    FC_GLOBAL_(slacpy,SLACPY)("All", &ctx->M, &ctx->M, ctx->RHS_X, &ctx->LDRHS_X, ctx->X, &ctx->LDX, 1);
#else
    FC_GLOBAL_(dlacpy,DLACPY)("All", &ctx->M, &ctx->M, ctx->RHS_X, &ctx->LDRHS_X, ctx->X, &ctx->LDX, 1);
#endif
}

static void context_lyap_iter_solve(void *_ctx){
    context_lyap_t * ctx = _ctx;
    int info= 0;

    if ( ctx->isolver < 7 ) {
        MEPACK_PREFIX(trlyap_level3)(ctx->TRANSA, ctx->M,
                ctx->A, ctx->LDA, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver < 12 ) {
        MEPACK_PREFIX(trlyap_dag)(ctx->TRANSA, ctx->M,
                ctx->A, ctx->LDA, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver < 19 ) {
        MEPACK_PREFIX(trlyap_level3_2stage)(ctx->TRANSA, ctx->M,
                ctx->A, ctx->LDA, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 19 ){
            MEPACK_PREFIX(trlyap_level2_opt)(ctx->TRANSA,  ctx->M,
                ctx->A, ctx->LDA, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 20 ) {
            MEPACK_PREFIX(trlyap_level2)(ctx->TRANSA, ctx->M,
                ctx->A, ctx->LDA, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    } else if ( ctx->isolver == 24){
        MEPACK_PREFIX(trlyap_recursive)(ctx->TRANSA, ctx->M,
                ctx->A, ctx->LDA, ctx->X, ctx->LDX,  &ctx->scale, ctx->work, &info);
    }
#if defined(RECSY) && !defined (SINGLE_PRECISION)
    else if ( ctx->isolver == 25 ) {
        FC_GLOBAL_(reclyct,RECLYCT)(&ctx->type, &ctx->scale, &ctx->M, ctx->A, &ctx->LDA, ctx->X, &ctx->LDX,  &info, ctx->MACHINE_RECSY);
    } else if ( ctx->isolver == 26 ) {
        Int proc = omp_get_num_procs();
        FC_GLOBAL_(reclyct_p,RECLYCT_P)(&proc, &ctx->type, &ctx->scale, &ctx->M, ctx->A, &ctx->LDA, ctx->X, &ctx->LDX, &info, ctx->MACHINE_RECSY);
    }
#endif
#ifndef SINGLE_PRECISION
    else if ( ctx->isolver == 30) {
        char TRANSX[2]={0,0};
        if ( tolower(ctx->TRANSA[0]) == 'n') {
            TRANSX[0]='T';
        } else {
            TRANSX[0]='N';
        }


        FC_GLOBAL_(sb03my,SB03MY)(TRANSX, &ctx->M, ctx->A, &ctx->LDA, ctx->X, &ctx->LDX, &ctx->scale, &info, 1);
    }
#endif
}

static void context_lyap_check(void *_ctx, double *forward, double *rel)
{
    context_lyap_t * ctx = _ctx;
    *rel += MEPACK_PREFIX(residual_lyap)(ctx->TRANSA, ctx->M,
            ctx->A, ctx->LDA,
            ctx->X, ctx->LDX, ctx->RHS_X, ctx->LDRHS_X, ctx->scale);
#ifdef SINGLE_PRECISION
    *forward += benchmark_check_X_float(ctx->M,ctx->M,ctx->X, ctx->LDX, ctx->Xorig, ctx->LDXorig);
#else
    *forward += benchmark_check_X_double(ctx->M,ctx->M,ctx->X, ctx->LDX, ctx->Xorig, ctx->LDXorig);
#endif
}


static void context_lyap_clean(void *_ctx)
{
    context_lyap_t * ctx = _ctx;

    free(ctx->A);
    free(ctx->X);
    free(ctx->Xorig);
    free(ctx->RHS_X);
    free(ctx->work);
    free(ctx->TRANSA);
}

benchmark_context_t *benchmark_create(void)
{
    benchmark_context_t *bctx = malloc(sizeof(benchmark_context_t));
    context_lyap_t * ctx = malloc(sizeof(context_lyap_t));
    memset(ctx, 0, sizeof(context_lyap_t));
    bctx->ctx = ctx;
    bctx->init = context_lyap_init;
    bctx->setup_rhs = context_lyap_setup_rhs;
    bctx->setup_problem = context_lyap_setup_problem;
    bctx->iter_prepare = context_lyap_iter_prepare;
    bctx->iter_solve = context_lyap_iter_solve;
    bctx->check = context_lyap_check;
    bctx->clean = context_lyap_clean;
    return bctx;
}

const char *equation_name(void) {
    return "Triangular Lyapunov Equation";
}

int have_solver(int is){
    if ( is > 0 && is <=26){
#ifndef RECSY
        if(is == 25 || is == 26 ) return 0;
#endif
#ifdef SINGLE_PRECISION
        if(is == 25 || is == 26 ) return 0;
#endif
        if (is == 21 || is == 22 || is == 23) return 0;
        return 1;
    }
#ifndef SINGLE_PRECISION
    if ( is == 30) return 1;
#endif
    return 0;

}
