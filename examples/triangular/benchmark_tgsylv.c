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

#ifdef RECSY
#define MAX_SOLVER 26
#else
#define MAX_SOLVER 24
#endif


void solver_name(int is) {
    switch(is) {
        case 1:
            printf("LEVEL3 - LEVEL2: LOCAL COPY fixed maximum size with alignment\n");
            break;
        case 2:
            printf("LEVEL3 - LEVEL2: LOCAL COPY fixed maximum size\n");
            break;
        case 3:
            printf("LEVEL3 - LEVEL2: REORDERED (column first solution)\n");
            break;
        case 4:
            printf("LEVEL3 - LEVEL2: BLAS LEVEL-2 call\n");
            break;
        case 5:
            printf("LEVEL3 - LEVEL2: Unoptimized\n");
            break;
        case 6:
            printf("LEVEL3 - LEVEL2: RECURSIVE BLOCKING\n");
            break;
        case 7:
            printf("DAG - LEVEL2: LOCAL COPY fixed maximum size with alignment\n");
            break;
        case 8:
            printf("DAG - LEVEL2: LOCAL COPY fixed maximum size\n");
            break;
        case 9:
            printf("DAG - LEVEL2: REORDERED (column first solution)\n");
            break;
        case 10:
            printf("DAG - LEVEL2: BLAS LEVEL-2 call\n");
            break;
        case 11:
            printf("DAG - LEVEL2: Unoptimized\n");
            break;
        case 12:
            printf("DAG - LEVEL2: RECURSIVE BLOCKING\n");
            break;
        case 13:
            printf("LEVEL3 - DAG - LEVEL2: LOCAL COPY fixed maximum size with alignment\n");
            break;
        case 14:
            printf("LEVEL3 - DAG - LEVEL2: LOCAL COPY fixed maximum size\n");
            break;
        case 15:
            printf("LEVEL3 - DAG - LEVEL2: REORDERED (column first solution)\n");
            break;
        case 16:
            printf("LEVEL3 - DAG - LEVEL2: BLAS LEVEL-2 call\n");
            break;
        case 17:
            printf("LEVEL3 - DAG - LEVEL2: Unoptimized\n");
            break;
        case 18:
            printf("LEVEL3 - DAG - LEVEL2: RECURSIVE BLOCKING\n");
            break;
        case 19:
            printf("LEVEL2: LOCAL COPY fixed maximum size with alignment\n");
            break;
        case 20:
            printf("LEVEL2: LOCAL COPY fixed maximum size\n");
            break;
        case 21:
            printf("LEVEL2: REORDERED (column first solution)\n");
            break;
        case 22:
            printf("LEVEL2: BLAS LEVEL-2 call\n");
            break;
        case 23:
            printf("LEVEL2: Unoptimized\n");
            break;
        case 24:
            printf("RECURSIVE BLOCKING\n");
            break;
        case 25:
            printf("RECSY\n");
            break;
        case 26:
            printf("RECSY Parallel\n");
            break;
        case 27:
            printf("Gardiner-Laub\n");
            break;
        case 28:
            printf("Gardiner-Laub (Algorithm 705)\n");
            break;

    }

}

void usage(char *prgmname) {
    int is;
    printf("Solve a generalized Sylvester equation with triangular coefficient matrices.\n");
    printf("\n");
    printf("Usage: %s <options>\n", prgmname);
    printf("\n");
    printf("--help, -h        Display this help\n");
    printf("--rows=M , -m M   Set the number of rows to M. \n");
    printf("--cols=N , -n N   Set the number of cols to N. \n");
    printf("--rows=START:STEP:STOP   Iterate the number of rows from START to STOP with step size STEP.\n");
    printf("--cols=START:STEP:STOP   Iterate the number of columns from START to STOP with step size STEP.\n");
    printf("--cols=-                 Use the number of rows as number of cols.\n");
    printf("--mb=START:STEP:STOP     Iterate the row block size form START to STOP with step size STEP.\n");
    printf("--nb=START:STEP:STOP     Iterate the column column block size from START to STOP with step size STEP.\n");
    printf("--mb=MB                  Set the row block size.  \n");
    printf("--nb=NB                  Set the column block size\n");
    printf("--nb=-                   Set the column block size to the row block size.\n");
    printf("--transa=N,T             Transpose on the matrix pair (A,C) \n");
    printf("--transb=N,T             Transpose on the matrix pair (B,D) \n");
    printf("--nmat=NUM               Number of different equations.\n");
    printf("--runs=R                 Number of runs per equation.\n");
    printf("--solver=S, -s S         Select the solver \n");
    printf("--alignoff, -a           Turn off the blocksize alignment.\n");
    printf("--changerole, -c         B and D changes places.\n");
    printf("--sign=+/-1, -S  +/-1    Sign in the equation.\n");
    printf("\nPossible Solvers\n");

    for (is = 0; is <= 28; is++) {
        printf("%2d : ", is); solver_name(is);
    }
    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i, mat ;
    double *A, *B, *C, *D;
    double *X, *Xorig,*Work, *RHS;
    Int M = 1024, N=1024;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int N_MIN=1024,N_MAX=1024,N_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int NB_MIN=64,NB_MAX=64,NB_STEP=32;
    Int MB, NB;

    Int RUNS = 5;
    Int is = 0;
    Int mn_same = 0;
    Int mbnb_same = 0;
    Int nMAT = 1;
    Int align_on = 1;
    double alpha, beta;
    double sign = 1;
    char TRANSA[20]="N";
    char TRANSB[20]="N";
#ifdef RECSY
    double MACHINE_RECSY[12];
    Int type, N2;
#endif

    double scale = 1.0;
    int info, run;
    double te, ts;
    ssize_t mem;

    double times,ts2, te2;
    double ctimes;
    double ress;
    double eps;
    int choice;
    int changerole = 0;

    NB = MB = 64;

    while (1)
    {
        static struct option long_options[] =
        {
            /* Use flags like so:
               {"verbose",    no_argument,    &verbose_flag, 'V'}*/
            /* Argument styles: no_argument, required_argument, optional_argument */
            {"help",    no_argument,    0,    'h'},
            {"rows",    required_argument, 0, 'm'},
            {"cols",    required_argument, 0, 'n'},
            {"mb",      required_argument, 0, 'M'},
            {"nb",      required_argument, 0, 'N'},
            {"transa",  required_argument, 0, 'A'},
            {"transb",  required_argument, 0, 'B'},
            {"nmat",    required_argument, 0, 't'},
            {"runs",    required_argument, 0, 'r'},
            {"solver",  required_argument, 0, 's'},
            {"changerole", no_argument, 0, 'c'},
            {"sign",  required_argument, 0, 'S'},
            {"alignoff", no_argument, 0, 'a'},

            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
no_argument: " "
required_argument: ":"
optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:n:M:N:A:B:t:r:s:S:ca",
                long_options, &option_index);

        if (choice == -1)
            break;

        switch( choice )
        {
            case 'h':
                usage(argv[0]);
                exit(-1);
                break;
            case 'm':
                if ( strstr(optarg, ":") != NULL ) {
                    int f1,f2,f3;
                    if ( sscanf(optarg, "%d:%d:%d", &f1, &f2, &f3) != 3 ) {
                        fprintf(stderr, "The loop argument must have the format MIN:STEP:MAX\n");
                    }
                    M_MIN = f1;
                    M_STEP = f2;
                    M_MAX = f3;
                    M = f1;
                } else {
                    M = atoi(optarg);
                    M_MIN = M;
                    M_MAX = M;
                    M_STEP = M;
                }
                break;
            case 'n':
                if ( strcmp("-", optarg) == 0 ) {
                    mn_same = 1;
                    N = M;
                    N_MIN = N;
                    N_MAX = N;
                    N_STEP = N;
                    break;
                }
                if ( strstr(optarg, ":") != NULL ) {
                    int f1,f2,f3;
                    if ( sscanf(optarg, "%d:%d:%d", &f1, &f2, &f3) != 3 ) {
                        fprintf(stderr, "The loop argument must have the format MIN:STEP:MAX\n");
                    }
                    N_MIN = f1;
                    N_STEP = f2;
                    N_MAX = f3;
                    N = f1;


                } else {
                    N = atoi(optarg);
                    N_MIN = N;
                    N_MAX = N;
                    N_STEP = N;
                }
                break;
            case 'M':
                if ( strstr(optarg, ":") != NULL ) {
                    int f1,f2,f3;
                    if ( sscanf(optarg, "%d:%d:%d", &f1, &f2, &f3) != 3 ) {
                        fprintf(stderr, "The loop argument must have the format MIN:STEP:MAX\n");
                    }
                    MB = f1;
                    MB_MIN = f1;
                    MB_STEP = f2;
                    MB_MAX = f3;
                } else {
                    MB_MIN = atoi(optarg);
                    MB_MAX = MB_MIN;
                    MB_STEP = MB_MIN;
                    MB = MB_MIN;
                }
                break;
            case 'N':
                if ( strcmp("-", optarg) == 0 ) {
                    mbnb_same = 1;
                    NB = MB;
                    NB_MIN = NB;
                    NB_MAX = NB;
                    NB_STEP = NB;
                    break;
                }
                if ( strstr(optarg, ":") != NULL ) {
                    int f1,f2,f3;
                    if ( sscanf(optarg, "%d:%d:%d", &f1, &f2, &f3) != 3 ) {
                        fprintf(stderr, "The loop argument must have the format MIN:STEP:MAX\n");
                    }
                    NB_MIN = f1;
                    NB_STEP = f2;
                    NB_MAX = f3;
                    NB = f1;


                } else {
                    NB = atoi(optarg);
                    NB_MIN = NB;
                    NB_MAX = NB;
                    NB_STEP = NB;
                }
                break;
            case 'A':
                strncpy(TRANSA, optarg,2);
                TRANSA[0] = toupper(TRANSA[0]);
                if (TRANSA[0] != 'N' && TRANSA[0] != 'T') {
                    printf("Invalid transpose option for (A,C).\n");
                    return -1;
                }
                break;
            case 'B':
                strncpy(TRANSB, optarg,2);
                TRANSB[0] = toupper(TRANSB[0]);
                if (TRANSB[0] != 'N' && TRANSB[0] != 'T') {
                    printf("Invalid transpose option for (B,D).\n");
                    return -1;
                }
                break;
            case 't':
                nMAT = atoi(optarg);
                break;
            case 'r':
                RUNS = atoi(optarg);
                break;
            case 's':
                is = atoi(optarg);
                break;
            case 'S':
                sign = atof(optarg);
                break;
            case 'a':
                align_on = 0;
                break;
            case 'c':
                changerole = 1;
                break;
            default:
                /* Not sure how to get here... */
                return EXIT_FAILURE;
        }
    }

    if ( is < 1 || is > 28 ) {
        fprintf(stderr, "Solver not known. \n");
        exit(-1);
    }

    benchmark_init();

    mepack_init();

    /*-----------------------------------------------------------------------------
     *  Output Configuration
     *-----------------------------------------------------------------------------*/
    printf("# Command Line: ");
    for (i = 1; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("# RUNS:  %d\n", (int) RUNS);
    printf("# Number of Matrices: %d\n",(int)  nMAT);
    printf("# Rows: %d (%d:%d:%d)\n", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    if ( mn_same )
        printf("# Cols: same as rows.\n");
    else
        printf("# Cols: %d (%d:%d:%d)\n", (int) N , (int) N_MIN, (int) N_STEP, (int) N_MAX);
    printf("# TRANSA: %s\n", TRANSA);
    printf("# TRANSB: %s\n", TRANSB);
    printf("# SIGN:   %lg\n", sign);
    printf("# Block Alignment: %s\n", (align_on == 1)?"YES":"NO");
    if ( changerole) printf("# B and D changes places.\n");
    mepack_tgsylv_isolver_set(1);
    mepack_double_tgsylv_blocksize_2stage_set(256);

    eps = mepack_double_epsilon();

    printf("# Solver: "); solver_name(is);


    printf("#\n");
    printf("#  M    N    MB  NB   Wall-Time     CPU-Time       Ratio    Forward-Err\n");
    for (M = M_MIN, N = N_MIN; M <= M_MAX && N<=N_MAX; M = (N>=N_MAX) ? (M+M_STEP): M,  N=(N<N_MAX)?(N+N_STEP):(N_MIN)) {
        if (mn_same) {
            N = M;
        }
        for (MB = MB_MIN, NB = NB_MIN; MB <= MB_MAX && NB<=NB_MAX; MB = (NB>=NB_MAX) ? (MB+MB_STEP): MB,  NB=(NB<NB_MAX)?(NB+NB_STEP):(NB_MIN)) {
            if ( mbnb_same ) {
                NB = MB;
            }

            /* Set the Block size   */
            mepack_double_tgsylv_blocksize_mb_set(MB);
            mepack_double_tgsylv_blocksize_nb_set(NB);

            /* Prepare  */
            times = 0;
            ctimes = 0;
            ress = 0.0;

            A = (double *) malloc(sizeof(double) * (M*M));
            C = (double *) malloc(sizeof(double) * (M*M));
            B = (double *) malloc(sizeof(double) * (N*N));
            D = (double *) malloc(sizeof(double) * (N*N));
            X = (double *) malloc(sizeof(double) * (M*N));
            memset(X, 0, sizeof(double) * M * N);
            /* #pragma omp parallel for schedule(static,1) */
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; C[i] = 0.0;  }
            /* #pragma omp parallel for schedule(static,1) */
            for (i = 0; i < N*N; i+=512) { B[i] = 0.0; D[i] = 0.0; }
            /* #pragma omp parallel for schedule(static,1) */
            for (i = 0; i < M*N; i+=1) { X[i] = 0.0; }

            Xorig = (double *) malloc(sizeof(double) * (M*N));
            RHS = (double *) malloc(sizeof(double) * (M*N));

            if ( is < 7 ) {
                 mepack_tgsylv_isolver_set(is);
                 mem = mepack_memory("DLA_TGSYLV_L3", M, N);
            } else if ( is > 6 && is < 13 ) {
                mepack_tgsylv_isolver_set(is-6);
                mem = mepack_memory("DLA_TGSYLV_DAG", M, N);
            } else if ( is > 12 && is < 19 ) {
                mepack_tgsylv_isolver_set(is-12);
                mem = mepack_memory("DLA_TGSYLV_L3_2S", M, N);
            } else if ( is == 19 ) {
                if ( M <= 32 && N <= 32)
                    mem = mepack_memory("DLA_TGSYLV_L2_LOCAL_COPY_32", M, N);
                else if ( M <= 64 && N <= 64)
                    mem = mepack_memory("DLA_TGSYLV_L2_LOCAL_COPY_64", M, N);
                else if ( M <= 96 && N <= 96)
                    mem = mepack_memory("DLA_TGSYLV_L2_LOCAL_COPY_96", M, N);
                else if ( M <= 128 && N <= 128)
                    mem = mepack_memory("DLA_TGSYLV_L2_LOCAL_COPY_128", M, N);
                else
                    mem = mepack_memory("DLA_TGSYLV_L2_REORDER", M, N);
            } else if ( is == 20 ) {
                if ( M <= 128 && N <= 128)
                    mem = mepack_memory("DLA_TGSYLV_L2_LOCAL_COPY", M, N);
                else
                    mem = mepack_memory("DLA_TGSYLV_L2_REORDER", M, N);
            }
            else if ( is == 20 ) {
                mem = mepack_memory("DLA_TGSYLV_L2_LOCAL_COPY", M, N);
            } else if ( is == 21) {
                mem = mepack_memory("DLA_TGSYLV_L2_REORDER", M, N);
            } else if ( is == 22) {
                mem = mepack_memory("DLA_TGSYLV_L2", M, N);
            } else if ( is == 23) {
                mem = mepack_memory("DLA_TGSYLV_L2_UNOPT", M, N);
            } else if ( is == 24) {
                mem = mepack_memory("DLA_TGSYLV_RECURSIVE", M, N);
            } else if ( is == 27) {
                mem = mepack_memory("DLA_TGSYLV_GARDINER_LAUB", M, N);
            } else {
                mem = 2*M*N;
            }

            if ( is == 19 || is == 20 ) {
                mem = MAX(mepack_memory("DLA_TGSYLV_L2_REORDER", M, N), mem);
            }
            if ( is == 28 ) {
                mem = 2*M*M+7*M;
            }
            Work = (double *) malloc(sizeof(double) * (mem));

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &N, &alpha, &beta, Xorig, &M);

            /* Reset the Seed  */
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                //		printf("seed %d %d %d %d\n", iseed[0], iseed[1], iseed[2], iseed[3]);
                benchmark_random_gevp_double(M, iseed, A, C, NULL, NULL, NULL, NULL, align_on*gcd(MB, NB));
                benchmark_random_gevp_double(N, iseed, B, D, NULL, NULL, NULL, NULL, align_on*gcd(MB, NB));
                // 		printf("seed %d %d %d %d\n", iseed[0], iseed[1], iseed[2], iseed[3]);

                if ( changerole ) {
                    double *tmp = malloc ( sizeof( double ) * N * N);
                    FC_GLOBAL(dlacpy,DLACPY) ( "All", &N, &N, B, &N, tmp, &N);
                    FC_GLOBAL(dlacpy,DLACPY) ( "All", &N, &N, D, &N, B, &N);
                    FC_GLOBAL(dlacpy,DLACPY) ( "All", &N, &N, tmp, &N, D, &N);
                    free(tmp);
                }
                benchmark_rhs_gsylv_double(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, Xorig, M, RHS, M );


                /*-----------------------------------------------------------------------------
                 *  LEVEL 3
                 *-----------------------------------------------------------------------------*/
                if ( is < 19 ) {
                    if ( is < 7 )  {
                        mepack_tgsylv_isolver_set(is);
                    } else if ( is > 6 && is < 13 )  {
                        mepack_tgsylv_isolver_set(is-6);
                    } else if ( is > 12 && is < 19 )  {
                        mepack_tgsylv_isolver_set(is-12);
                    }

                    te = 0.0;
                    te2 = 0.0;
                    for (run = -1; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M);
                        ts = get_wtime();
                        ts2 = get_ctime();
                        if ( is < 7 )
                            mepack_double_tgsylv_level3(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        else if ( is < 12 )
                            mepack_double_tgsylv_dag(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        else if ( is < 19 )
                            mepack_double_tgsylv_level3_2stage(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);

                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }

                    }
                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;

                    ress += benchmark_check_X_double(M,N,X, M, Xorig, M);
                }

                /*-----------------------------------------------------------------------------
                 *  LEVEL 2
                 *-----------------------------------------------------------------------------*/
                if ( is > 18 && is < 25 )  {
                    te = 0.0;
                    te2 = 0.0;


                    for (run = -1; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M);
                        ts = get_wtime();
                        ts2 = get_ctime();
                        if ( is == 19 ){
                            if ( M <= 32 && N <= 32 )
                                mepack_double_tgsylv_level2_local_copy_32(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                            else if ( M <= 64 && N <= 64 )
                                mepack_double_tgsylv_level2_local_copy_64(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                            else if ( M <= 96 && N <= 96 )
                                mepack_double_tgsylv_level2_local_copy_96(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                            else if ( M <= 128 && N <= 128 )
                                mepack_double_tgsylv_level2_local_copy_128(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                            else
                                mepack_double_tgsylv_level2_reorder(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        } else if ( is == 20 ) {
                            if ( M <= 128 && N <= 128) {
                                mepack_double_tgsylv_level2_local_copy(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                            } else {
                                mepack_double_tgsylv_level2_reorder(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                            }
                        } else if ( is == 21 ) {
                            mepack_double_tgsylv_level2_reorder(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        } else if ( is == 22 ) {
                            mepack_double_tgsylv_level2(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        } else if ( is == 23 ) {
                            mepack_double_tgsylv_level2_rowwise(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        } else {
                            mepack_double_tgsylv_recursive(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        }
                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }
                    }


                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;


                    ress += benchmark_check_X_double(M,N,X, M, Xorig, M);
                }
#ifdef RECSY
                /*-----------------------------------------------------------------------------
                 *  RECSY
                 *-----------------------------------------------------------------------------*/
                type = benchmark_sylv_setup_to_type(TRANSA, TRANSB, sign);
                FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(MACHINE_RECSY);
                N2 = M*N;
                if ( is == 25  && changerole == 0)  {
                    Int infox= 0;
                    te = 0.0;
                    te2 = 0.0;
                    for (run = -1; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M);
                        ts = get_wtime();
                        ts2 = get_ctime();
                        FC_GLOBAL_(recgsyl,RECGSYL)(&type, &scale, &M, &N, A, &M, B, &N,  C, &M, D, &N, X,  &M, &infox, MACHINE_RECSY, Work, &N2);

                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }
                    }
                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;


                    ress += benchmark_check_X_double(M,N,X, M, Xorig, M);
                }


                /*-----------------------------------------------------------------------------
                 *  RECSY_P
                 *-----------------------------------------------------------------------------*/
                if ( is == 26 && changerole == 0 )  {
                    Int infox = 0;
                    Int proc = omp_get_max_threads();
                    te = 0.0;
                    te2 = 0.0;
                    for (run = -1; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M);
                        ts = get_wtime();
                        ts2 = get_ctime();
                        FC_GLOBAL_(recgsyl_p,RECGSYL_P)(&proc,&type, &scale, &M, &N, A, &M, B, &N,  C, &M, D, &N, X,  &M, &infox, MACHINE_RECSY, Work, &N2);

                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }
                    }
                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;


                    ress += benchmark_check_X_double(M,N,X, M, Xorig, M);
                }


#endif
                /*-----------------------------------------------------------------------------
                 *  Gardiner Laub
                 *-----------------------------------------------------------------------------*/
                if ( is == 27 )  {
                    te = 0.0;
                    te2 = 0.0;
                    for (run = -1; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M);
                        ts = get_wtime();
                        ts2 = get_ctime();
                        mepack_double_tgsylv_gardiner_laub(TRANSA, TRANSB, sign,  M, N, A, M, B, N, C, M, D, N, X, M, &scale, Work, &info);
                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }

                    }
                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;

                    ress += benchmark_check_X_double(M,N,X, M, Xorig, M);
                }

#ifdef GL705
                /*-----------------------------------------------------------------------------
                 *  Gardiner Laub (Old Algorithm 705 Code)
                 *-----------------------------------------------------------------------------*/
                if ( is == 28 )  {
                    Int job = 0;
                    Int *IWork;
                    if (tolower(TRANSA[0])!='n' && tolower(TRANSB[0]) != 't') {
                        printf("# Gardiner-Laub (Algorithm 705) only supports the (N,T) case\n");
                    }
                    if (!changerole) {
                        printf("# Gardiner-Laub (Algorithm 705) only works if B and D change roles.\n");
                    }
                    IWork = (Int *) malloc(sizeof(Int) * (2*M));
                    Int infox = 0;
                    te = 0.0;
                    te2 = 0.0;
                    for (run = -1; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M);
                        ts = get_wtime();
                        ts2 = get_ctime();
                        // FC_GLOBAL_(dla_tgsylv_gardiner_laub,DLA_TGSYLV_GARDINER_LAUB)(TRANSA, TRANSB, &sign,  &M, &N, A, &M, B, &N, C, &M, D, &N,  X, &M, &scale, Work, &info, MACHINE);
                        FC_GLOBAL(bkhs2,BKHS2) (&M, &N, A, &M, B, &N, C,&M,  D,&N,  X,&M,  Work, IWork, &job, &infox);
                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }

                    }
                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;
                    free(IWork);
                    ress += benchmark_check_X_double(M,N,X, M, Xorig, M);
                }
#endif
            }  // nmat
            times /= (double) nMAT;
            ctimes /= (double) nMAT;
            ress /= (double) nMAT;

            /* Print  */
            printf("%5d %5d %3d %3d  %10.5e  %10.5e  %10.5e  %10.5e \n",
                   (int)  M,(int) N, (int) MB, (int) NB, times, ctimes, ctimes/times, ress);
            fflush(stdout);
            free(A);
            free(B);
            free(C);
            free(D);
            free(X);
            free(Xorig);
            free(Work);
            free(RHS);

        }
    }

    eps = mepack_double_epsilon();

    benchmark_exit();
    /* printf("ress %lg %lg\n", ress, sqrt(eps)*100); */
    return (ress < sqrt(eps)*100)? 0 : -1 ;

}


