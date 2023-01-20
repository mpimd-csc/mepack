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
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <sys/time.h>
#include <stdint.h>
#include <time.h>
#include <complex.h>
#include "cscutils/io.h"
#include "cscutils/hdf.h"

#include "benchmark.h"


char *hdfstore_path = NULL;

void benchmark_init(void)
{
    char * tmp = malloc(sizeof(char) * PATH_LEN + 128);
    hdfstore_path = malloc(sizeof(char) * PATH_LEN + 128);

    if (getenv("MEPACK_HDF_PATH") == NULL) {
        strncpy(hdfstore_path, "./", PATH_LEN);
    } else {
        strncpy(hdfstore_path, getenv("MEPACK_HDF_PATH"), PATH_LEN);
    }
    fprintf(stdout, "# HDF5 Store Path: %s (set with MEPACK_HDF_PATH)\n", hdfstore_path);

    snprintf(tmp, PATH_LEN, "%s/qr", hdfstore_path);
    if ( csc_file_mkdir(tmp, S_IRUSR|S_IWUSR|S_IXUSR) ) {
        fprintf(stderr, "Failed to create %s\n", tmp);
        exit(-1);
    }

    snprintf(tmp, PATH_LEN, "%s/qz", hdfstore_path);
    if ( csc_file_mkdir(tmp, S_IRUSR|S_IWUSR|S_IXUSR)) {
        fprintf(stderr, "Failed to create %s\n", tmp);
        exit(-1);
    }
    free(tmp);

    return ;

}

void FC_GLOBAL(benchmark_init_f,BENCHMARK_INIT_F) (void)
{
    benchmark_init();
}

void benchmark_exit(void)
{
    free(hdfstore_path);
    return;
}

double get_wtime(void)
{
	struct timeval tv;
	gettimeofday (&tv, NULL);
	return tv.tv_sec + tv.tv_usec / 1e6;
}

double get_ctime(void) {
	clock_t  t = clock();
	return t / (double)(CLOCKS_PER_SEC);
}

Int gcd(Int a, Int b)
{
    Int h;
    while ( b != 0 ) {
        h = a % b;
        a = b;
        b = h;
    }
    return a;
}


void benchmark_print_matrix_double(Int M, Int N, double *X, Int LDX)
{
    Int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            printf("%6.5e \t", X[i+j*LDX]);
        }
        printf("\n");
    }
}

void benchmark_print_matrix_float(Int M, Int N, float *X, Int LDX)
{
    Int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            printf("%6.5e \t", (double) X[i+j*LDX]);
        }
        printf("\n");
    }
}


__attribute__((weak))
void solver_name(int is) {
     switch(is) {
        case MEPACK_FRONTEND_SOLVER_LEVEL3:
            printf("Refine with LEVEL3 solver\n");
            break;
        case MEPACK_FRONTEND_SOLVER_LEVEL2:
            printf("Refine with LEVEL2 solver\n");
            break;
        case MEPACK_FRONTEND_SOLVER_2STAGE:
            printf("Refine with LEVEL3 two stage solver\n");
            break;
        case MEPACK_FRONTEND_SOLVER_DAG:
            printf("Refine with DAG solver.\n");
            break;
        case MEPACK_FRONTEND_SOLVER_RECURSIVE:
            printf("Refine with RECURSIVE solver.\n");
            break;
    }
}

