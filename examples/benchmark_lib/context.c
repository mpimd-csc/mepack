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

void benchmark_context_free(benchmark_context_t *bctx)
{
    free(bctx->ctx);
    free(bctx);
}

void benchmark_context_init(benchmark_context_t *ctx, Int M, Int N, Int ISolver, Int MB, Int NB, Int BIGMB, const char *TA, const char *TB, double sgn1, double sgn2)
{
    ctx->init(ctx->ctx, M, N, ISolver, MB, NB, BIGMB, TA, TB, sgn1, sgn2);
}

void benchmark_context_setup_rhs(benchmark_context_t  *ctx){
    ctx->setup_rhs(ctx->ctx);
}

void benchmark_context_setup_problem(benchmark_context_t *ctx, Int iseed[4], Int align, Int changerole)
{
    ctx->setup_problem(ctx->ctx, iseed, align, changerole);

}

void benchmark_context_iter_prepare(benchmark_context_t *ctx)
{
    ctx->iter_prepare(ctx->ctx);

}

void benchmark_context_iter_solve(benchmark_context_t *ctx){
    ctx->iter_solve(ctx->ctx);
}

void benchmark_context_check(benchmark_context_t *ctx, double *forward, double *rel)
{
    ctx->check(ctx->ctx, forward, rel);
}

void benchmark_context_clean(benchmark_context_t *ctx)
{
    ctx->clean(ctx->ctx);

}

