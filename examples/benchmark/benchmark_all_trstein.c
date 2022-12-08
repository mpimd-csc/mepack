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

#include "benchmark.h"
#include "mepack.h"
#include "mepack_internal.h"



void usage(char *prgmname) {
    printf("Benchmark the different level 2 inside the level 3 solver for the standard Stein equation.\n");
    printf("\n");
    printf("Usage: %s <options>\n",prgmname);
    printf("\n");
    printf("--help, -h        Display this help\n");
    printf("--rows=M , -m M   Set the number of rows to M. \n");
    printf("--rows=START:STEP:STOP   Iterate the number of rows from START to STOP with step size STEP.\n");
    printf("--mb=START:STEP:STOP     Iterate the row block size form START to STOP with step size STEP.\n");
    printf("--mb=MB                  Set the row block size.  \n");
    printf("--trans=N,T              Transpose on the matrix A.\n");
    printf("--nmat=NUM               Number of different equations.\n");
    printf("--runs=R                 Number of runs per equation.\n");
    printf("--mb=START:STEP:STOP     Iterate over an interval of row blocks sizes.\n");
    printf("--mb=MB                  Set the row block size.  \n");
    printf("--res                    Compute and print the residual instead of the run time. \n");

    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i,mat ;
    double *A;
    double *X, *Xorig,*Work, *RHS;
    Int M = 1024;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int MB;

    Int RUNS = 5;
    Int is = 0;
    Int nMAT = 1;

    double alpha, beta;
    char TRANSA[20]="N";
    Int print_res = 0;
#ifdef RECSY
    double MACHINE_RECSY[12];
    Int type, j;
#endif

    double scale = 1.0;
    int info;
    Int run;
    double te, ts;


    double times[9];
    double ress[9];

    int choice;

    MB = 64;

    while (1)
    {
        static struct option long_options[] =
        {
            /* Use flags like so:
               {"verbose",    no_argument,    &verbose_flag, 'V'}*/
            /* Argument styles: no_argument, required_argument, optional_argument */
            {"help",    no_argument,    0,    'h'},
            {"rows",    required_argument, 0, 'm'},
            {"mb",      required_argument, 0, 'M'},
            {"trans",  required_argument, 0, 'A'},
            {"nmat",    required_argument, 0, 't'},
            {"runs",    required_argument, 0, 'r'},
            {"res",     no_argument, 0, 'R'},
            {0,0,0,0}
        };

        int option_index = 0;


        choice = getopt_long( argc, argv, "hm:M:A:t:r:R",
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
            case 'A':
                strncpy(TRANSA, optarg,2);
                TRANSA[0] = toupper(TRANSA[0]);
                if (TRANSA[0] != 'N' && TRANSA[0] != 'T') {
                    printf("Invalid transpose option for (A,C).\n");
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
    printf("# Rows: %d (%d:%d:%d)\n",(int)  M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    printf("# TRANSA: %s\n", TRANSA);
    if (print_res )  printf("# Print Residual\n");

    mepack_trstein_isolver_set(1);
   if ( getenv("MEPACK_LUA_CONFIG")) {
        printf("# MEPACK_LUA_CONFIG set. Using blocksize from %s\n", getenv("MEPACK_LUA_CONFIG"));
        MB_MIN = 1;
        MB_MAX = 1;
    }


    printf("#\n");
    printf("#  M   MB       RECSY     LC_ALIGN   LOCAL_COPY      REORDER       L2-OPT        UNOPT    RECURSIVE     OpenMP-DAG    L3-2Stage\n");
    for (M = M_MIN ;  M <= M_MAX ; M = M + M_STEP ) {
        for (MB = MB_MIN; MB <= MB_MAX ; MB += MB_STEP) {
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            /* Set the Block size   */
            if ( getenv("MEPACK_LUA_CONFIG")) {
                mepack_double_trstein_blocksize_set(0);
            }
            else
            {
                mepack_double_trstein_blocksize_set(MB);
            }



            for (i = 0; i < 9; i++) {
                times[i] = 0;
                ress[i] = 0.0;
            }


            A = (double *) malloc(sizeof(double) * (M*M));
            X = (double *) malloc(sizeof(double) * (M*M));
            Xorig = (double *) malloc(sizeof(double) * (M*M));
            RHS = (double *) malloc(sizeof(double) * (M*M));
            Work = (double *) malloc(sizeof(double) * (2*M*M));

#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0;}
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { X[i] = 0.0; }

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &M, &alpha, &beta, Xorig, &M);


            for (mat = 0; mat < nMAT; mat++) {
                benchmark_random_evp_double(M, iseed, A,  NULL, NULL, MB);
                benchmark_rhs_stein_double (TRANSA,  M, A, M, Xorig, M, RHS, M );

#ifdef RECSY
                /*-----------------------------------------------------------------------------
                 *  RECSY
                 *-----------------------------------------------------------------------------*/
                Int infox = 0;
                type = benchmark_gstein_setup_to_type(TRANSA);
                FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(MACHINE_RECSY);
                te = 0;
                j = M * M ;
                for (run = 0; run < RUNS; run++) {
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, RHS, &M, X, &M);
                    ts = get_wtime();
                    FC_GLOBAL_(reclydt,RECLYDT)(&type, &scale, &M, A, &M, X, &M, &infox, MACHINE_RECSY, Work, &j);
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < j; i++) {
                            X[i*M+j] = X[j*M+i];
                        }
                    }
                    te += (get_wtime() - ts);
                }
                te = te / RUNS;
                times[0] += te;
                ress[0] += benchmark_check_X_double(M,M,X,M,Xorig,M);
#else
                times[0] = 0;
                ress[0] = 0;
#endif
                /*-----------------------------------------------------------------------------
                 *  Bartels-Stewart
                 *-----------------------------------------------------------------------------*/
                for (is = 1; is < 9; is++) {
                    mepack_trstein_isolver_set(is);
                    te = 0.0;
                    for (run = 0; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, RHS, &M, X, &M);
                        if ( is < 7 ) {
                            ts = get_wtime();
                            mepack_double_trstein_level3(TRANSA,  M, A, M, X, M, &scale, Work, &info);
                            te += (get_wtime()-ts);
                        } else if ( is == 7 ) {
                            mepack_trstein_isolver_set(1);
                            ts = get_wtime();
                            mepack_double_trstein_dag(TRANSA,  M, A, M, X, M, &scale, Work, &info);
                            te += (get_wtime()-ts);
                        } else if ( is == 8 ) {
                            mepack_trstein_isolver_set(1);
                            ts = get_wtime();
                            mepack_double_trstein_level3_2stage(TRANSA,  M, A, M, X, M, &scale, Work, &info);
                            te += (get_wtime()-ts);
                        }
                    }
                    te = te / RUNS;
                    times[is] += te;
                    ress[is] += benchmark_check_X_double(M,M,X,M,Xorig,M);
                }
            }
            for (i = 0; i < 9; i++) {
                times[i] /= (double) nMAT;
                ress[i] /= (double) nMAT;
            }

            /* Print  */

            if (print_res ) {
                printf("%5d %3d  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e\n",
                        (int) M, (int) MB, ress[0], ress[1], ress[2], ress[3], ress[4], ress[5], ress[6] , ress[7], ress[8] );

            } else {
                printf("%5d %3d  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e  %10.5e\n",
                        (int) M,(int) MB, times[0], times[1], times[2], times[3], times[4], times[5], times[6], times[7], times[8] );
            }
            fflush(stdout);
            free(A);
            free(X);
            free(Xorig);
            free(Work);
            free(RHS);

        }
    }

    benchmark_exit();
    return 0;
}

