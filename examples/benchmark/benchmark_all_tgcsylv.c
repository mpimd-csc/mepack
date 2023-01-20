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

#include "benchmark.h"
#include "mepack.h"
#include "mepack_internal.h"

void usage(char * prgmname) {
    printf("Benchmark the different level 2 inside the level 3 solver for the generalized coupled Sylvester equation.\n");
    printf("\n");
    printf("Usage: %s <options>\n",prgmname);
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
    printf("--res                    Compute and print the residual instead of the run time. \n");
    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i,mat ;
    double *A, *B, *C, *D;
    double *Work;
    double *R, *Rorig, *RHS_R;
    double *L, *Lorig, *RHS_L;

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

    double alpha, beta;
    double sign1 = 1;
    double sign2 = 1;
    char TRANSA[20]="N";
    char TRANSB[20]="N";
    Int print_res = 0;
#ifdef RECSY
    double MACHINE_RECSY[12];
    Int type ;
#endif

    double scale = 1.0;
    int info, run;
    double te, ts;


    double times[9];
    double ress[9];

    int choice;

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
            {"res",     no_argument, 0, 'R'},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
no_argument: " "
required_argument: ":"
optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:n:M:N:A:B:t:r:R",
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
            case 'R':
                print_res =1;
                RUNS = 1;
                break;
            default:
                /* Not sure how to get here... */
                return EXIT_FAILURE;
        }
    }

    mepack_init();
    benchmark_init();
    /*-----------------------------------------------------------------------------
     *  Output Configuration
     *-----------------------------------------------------------------------------*/
    printf("# Command Line: ");
    for (i = 1; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("# RUNS:  %d\n", (int) RUNS);
    printf("# Number of Matrices: %d\n", (int) nMAT);
    printf("# Rows: %d (%d:%d:%d)\n", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    if ( mn_same )
        printf("# Cols: same as rows.\n");
    else
        printf("# Cols: %d (%d:%d:%d)\n", (int) N , (int) N_MIN, (int) N_STEP,(int)  N_MAX);
    printf("# TRANSA: %s\n", TRANSA);
    printf("# TRANSB: %s\n", TRANSB);
    if (print_res )  printf("# Print Residual\n");


    mepack_tgcsylv_isolver_set(1);

    if ( getenv("MEPACK_LUA_CONFIG")) {
        printf("# MEPACK_LUA_CONFIG set. Using blocksize from %s\n", getenv("MEPACK_LUA_CONFIG"));
        MB_MIN = 1;
        MB_MAX = 1;
        NB_MIN = 1;
        NB_MAX = 1;
    }


    printf("#\n");
    printf("#  M    N    MB  NB       RECSY     LC_ALIGN   LOCAL_COPY      REORDER       L2-OPT        UNOPT    RECURSIVE    OpenMP-DAG   L3-2STAGE\n");
    for (M = M_MIN, N = N_MIN; M <= M_MAX && N<=N_MAX; M = (N>=N_MAX) ? (M+M_STEP): M,  N=(N<N_MAX)?(N+N_STEP):(N_MIN)) {
        if (mn_same) {
            N = M;
        }
        for (MB = MB_MIN, NB = NB_MIN; MB <= MB_MAX && NB<=NB_MAX; MB = (NB>=NB_MAX) ? (MB+MB_STEP): MB,  NB=(NB<NB_MAX)?(NB+NB_STEP):(NB_MIN)) {
            if ( mbnb_same ) {
                NB = MB;
            }
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            /* Set the Block size   */

            if ( getenv("MEPACK_LUA_CONFIG")) {
                mepack_double_tgcsylv_blocksize_mb_set(0);
                mepack_double_tgcsylv_blocksize_nb_set(0);
            }
            else
            {
                mepack_double_tgcsylv_blocksize_mb_set(MB);
                mepack_double_tgcsylv_blocksize_nb_set(NB);
            }


            for (i = 0; i < 9; i++) {
                times[i] = 0;
                ress[i] = 0.0;
            }

            A = (double *) malloc(sizeof(double) * (M*M));
            C = (double *) malloc(sizeof(double) * (M*M));
            B = (double *) malloc(sizeof(double) * (N*N));
            D = (double *) malloc(sizeof(double) * (N*N));
            R = (double *) malloc(sizeof(double) * (M*N));
            L = (double *) malloc(sizeof(double) * (M*N));


#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; C[i] = 0.0;  }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < N*N; i+=512) { B[i] = 0.0; D[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*N; i+=512) { R[i] = 0.0; L[i] = 0.0;  }

            Rorig = (double *) malloc(sizeof(double) * (M*N));
            Lorig = (double *) malloc(sizeof(double) * (M*N));
            RHS_R = (double *) malloc(sizeof(double) * (M*N));
            RHS_L = (double *) malloc(sizeof(double) * (M*N));

            Work = (double *) malloc(sizeof(double) * (2*M*N));

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &N, &alpha, &beta, Rorig, &M, 1);
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &N, &alpha, &beta, Lorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                benchmark_random_gevp_double(M, iseed, A, C, NULL, NULL, NULL, NULL, gcd(MB, NB));
                benchmark_random_gevp_double(N, iseed, B, D, NULL, NULL, NULL, NULL, gcd(MB, NB));
                benchmark_rhs_ggcsylv_double(TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N, Rorig, M, Lorig, M,  RHS_R, M, RHS_L, M );

#ifdef RECSY
                /*-----------------------------------------------------------------------------
                 *  RECSY
                 *-----------------------------------------------------------------------------*/
                Int infox = 0;
                type = benchmark_ggcsylv_setup_to_type(TRANSA, TRANSB, sign1);
                FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(MACHINE_RECSY);
                te = 0;
                for (run = 0; run < RUNS; run++) {
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS_R, &M, R, &M, 1);
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS_L, &M, L, &M, 1);

                    ts = get_wtime();
                    FC_GLOBAL_(recgcsy,RECGCSY)(&type, &scale, &M, &N, A, &M, B, &N, R, &M, C, &M, D, &N, L, &M,  &infox, MACHINE_RECSY);
                    te += (get_wtime() - ts);
                }
                te = te / RUNS;
                times[0] += te;

                ress[0] += benchmark_check_X_double(M,N,R,M,Rorig,M);
                ress[0] += benchmark_check_X_double(M,N,L,M,Lorig,M);
#else
                times[0] = 0;
                ress[0] = 0;
#endif

                /*-----------------------------------------------------------------------------
                 *  Bartels-Stewart
                 *-----------------------------------------------------------------------------*/
                for (is = 1; is < 9; is++) {
                    mepack_tgcsylv_isolver_set(is);
                    te = 0.0;
                    for (run = 0; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS_R, &M, R, &M, 1);
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS_L, &M, L, &M, 1);

                        if ( is < 7 ) {
                            ts = get_wtime();
                            mepack_double_tgcsylv_level3(TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N, R, M, L, M,  &scale, Work, &info);
                            te += (get_wtime()-ts);
                        } else if ( is == 7 ) {
                            mepack_tgcsylv_isolver_set(1);
                            ts = get_wtime();
                            mepack_double_tgcsylv_dag(TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N, R, M, L, M,  &scale, Work, &info);
                            te += (get_wtime()-ts);
                        } else if ( is == 8 ) {
                            mepack_tgcsylv_isolver_set(1);
                            ts = get_wtime();
                            mepack_double_tgcsylv_level3_2stage(TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N, R, M, L, M,  &scale, Work, &info);
                            te += (get_wtime()-ts);
                        }
                    }
                    te = te / RUNS;
                    times[is] += te;
                    ress[is] += benchmark_check_X_double(M,N,R,M,Rorig,M);
                    ress[is] += benchmark_check_X_double(M,N,L,M,Lorig,M);
                }
            }
            for (i = 0; i < 9; i++) {
                times[i] /= (double) nMAT;
                ress[i] /= (double) nMAT;

            }

            /* Print  */

            if (print_res ) {
                printf("%5d %5d %3d %3d  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e %10.5e  %10.5e\n",
                        (int) M,(int) N, (int) MB, (int) NB, ress[0], ress[1], ress[2], ress[3], ress[4], ress[5], ress[6], ress[7], ress[8]);

            } else {
                printf("%5d %5d %3d %3d  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e\n",
                        (int) M,(int) N, (int) MB, (int) NB, times[0], times[1], times[2], times[3], times[4], times[5], times[6], times[7], times[8]);
            }
            fflush(stdout);
            free(A);
            free(B);
            free(C);
            free(D);
            free(R);
            free(Rorig);
            free(RHS_R);
            free(L);
            free(Lorig);
            free(RHS_L);
            free(Work);

        }
    }

        benchmark_exit();
    return 0;
}


