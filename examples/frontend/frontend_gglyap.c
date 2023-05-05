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
            printf("LEVEL2: Optimized DSYR2 Call\n");
            break;
        case 20:
            printf("Not Implemented\n");
            break;
        case 21:
            printf("Not Implemented\n");
            break;
        case 22:
            printf("Not Implemented\n");
            break;
        case 23:
            printf("Not Implemented\n");
            break;
        case 24:
            printf("RECURSIVE BLOCKING\n");
            break;
        case 27:
            printf("SLICOT\n");
            break;
            /* case 25:
               printf("RECSY\n");
               break;
               case 26:
               printf("RECSY Parallel\n");
               break;
               case 27:
               printf("SLICOT\n");
               break;
               case 28:
               printf("Gardiner-Laub (Algorithm 705)\n");
               break;  */
    }

}

void usage(char *prgmname) {
    int is;
    printf("Solve a generalized Lyapunov equation with different solvers.\n");
    printf("\n");
    printf("Usage: %s <options>\n",prgmname);
    printf("\n");
    printf("--help, -h        Display this help\n");
    printf("--rows=M , -m M   Set the number of rows to M. \n");
    printf("--rows=START:STEP:STOP   Iterate the number of rows from START to STOP with step size STEP.\n");
    printf("--mb=START:STEP:STOP     Iterate the row block size form START to STOP with step size STEP.\n");
    printf("--mb=MB                  Set the row block size.  \n");
    printf("--trans=N,T              Transpose on the matrix pair (A,C) \n");
    printf("--nmat=NUM               Number of different equations.\n");
    printf("--runs=R                 Number of runs per equation.\n");
    printf("--solver=S, -s S         Select the solver \n");
    printf("--alignoff, -a           Turn off the block alignment\n");
    printf("--reuse, -R              Benchmark the Reuse of the schur decomposition.\n");
    printf("\nPossible Solvers\n");
    for (is = 0; is <= 27; is++) {
        printf("%2d : ", is);
        solver_name(is);
    }

    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i,mat ;
    double *A, *B;
    double *Aorig, *Borig;
    double *Q, *Z;
    double *X, *Xorig,*Work, *RHS;
    Int M = 1024, N2;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int MB;

    Int RUNS = 5;
    Int is = 0;
    Int nMAT = 1;
    Int alignon = 1;
    size_t mem = 0 ;

    double alpha, beta;
    char TRANSA[20]="N";
    char reusestr[10];
    double eps;

    double scale = 1.0;
    int info, run;
    double te, ts = 0;

    double times,ts2 = 0, te2;
    double ctimes;
    double ress = 1.0;

    int choice;
    int reuse = 0;

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
            {"solver",  required_argument, 0, 's'},
            {"alignoff", no_argument, 0, 'a'},
            {"reuse", no_argument, 0, 'R',},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
no_argument: " "
required_argument: ":"
optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:M:A:t:r:s:S:R",
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
            case 's':
                is = atoi(optarg);
                break;
            case 'a':
                alignon = 0;
                break;
            case 'R':
                reuse = 1;
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
    printf("# Number of Matrices: %d\n",(int)  nMAT);
    printf("# Rows: %d (%d:%d:%d)\n", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    printf("# TRANSA: %s\n", TRANSA);
    printf("# Block Alignment: %s\n", (alignon == 1)?"YES":"NO");

    mepack_tglyap_isolver_set(1);
    mepack_double_tglyap_blocksize_2stage_set(256);

    eps = mepack_double_epsilon();

    printf("# Solver: "); solver_name(is);


    printf("#\n");
    printf("#  M   MB  Wall-Time     CPU-Time       Ratio    Forward-Err\n");

    for (M = M_MIN ;  M <= M_MAX ; M = M + M_STEP ) {
        for (MB = MB_MIN; MB <= MB_MAX ; MB += MB_STEP) {
            /* Reset the Seed  */
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;


            /* Set the Block size   */
            mepack_double_tglyap_blocksize_set(MB);

            /* Prepare  */
            times = 0;
            ctimes = 0;
            ress = 0.0;

            A = (double *) malloc(sizeof(double) * (M*M));
            B = (double *) malloc(sizeof(double) * (M*M));
            Q = (double *) malloc(sizeof(double) * (M*M));
            Z = (double *) malloc(sizeof(double) * (M*M));
            Aorig = (double *) malloc(sizeof(double) * (M*M));
            Borig = (double *) malloc(sizeof(double) * (M*M));
            X = (double *) malloc(sizeof(double) * (M*M));
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; B[i] = 0.0; Q[i] =0.0; Z[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { X[i] = 0.0; }

            Xorig = (double *) malloc(sizeof(double) * (M*M));
            RHS = (double *) malloc(sizeof(double) * (M*M));


            if ( is < 7 )  {
                mepack_tglyap_isolver_set(is);
                mepack_tglyap_frontend_solver_set(MEPACK_FRONTEND_SOLVER_LEVEL3);
            } else if ( is > 6 && is < 13 )  {
                mepack_tglyap_isolver_set(is-6);
                mepack_tglyap_frontend_solver_set(MEPACK_FRONTEND_SOLVER_DAG);
            } else if ( is > 12 && is < 19 )  {
                mepack_tglyap_isolver_set(is-12);
                mepack_tglyap_frontend_solver_set(MEPACK_FRONTEND_SOLVER_2STAGE);
            } else if ( is == 24  ) {
                mepack_tglyap_isolver_set(1);
                mepack_tglyap_frontend_solver_set(MEPACK_FRONTEND_SOLVER_RECURSIVE);
            }

            if ( reuse ) {
                strcpy(reusestr, "Nofact");
            } else {
                strcpy(reusestr, "Fact");
            }
            mem = mepack_memory_frontend("DLA_GGLYAP", reusestr, reusestr, M, M);



            Work = (double *) malloc(sizeof(double) * (mem));

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &M, &alpha, &beta, Xorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                // benchmark_random_gevp_double(M, iseed, A, B, NULL, NULL, NULL, NULL, alignon*MB);
                Int IDIST = 2;
                Int ldwork = mem;
                char TRANSX[2] = {0,0};
                double sep, ferr;
                Int ldworkx = ldwork - 3*M;
                N2 = M * M;
                FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, A);
                FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, B);
                FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, A, &M, Aorig, &M, 1);
                FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, B, &M, Borig, &M, 1);
                benchmark_rhs_glyap_double(TRANSA, M, A, M, B, M, Xorig, M, RHS, M );


                if ( tolower(TRANSA[0]) == 'n') {
                    TRANSX[0]='T';
                } else {
                    TRANSX[0]='N';
                }


                te = 0.0;
                te2 = 0.0;
                for (run = -1; run < RUNS; run++) {
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, RHS, &M, X, &M, 1);
                    if ( run == -1 || !reuse ) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, Aorig, &M, A, &M, 1);
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, Borig, &M, B, &M, 1);
                    }

                    ts = get_wtime();
                    ts2 = get_ctime();

                    if ( run == -1 || !reuse) {
                        if ( is == 27) {
                            Int infox = 0 ;
                            FC_GLOBAL_(sg03ad,SG03AD)("Cont", "X", "NonFact", TRANSX, "Upper", &M,
                                    A, &M, B, &M, Q, &M, Z, &M, X, &M, &scale, &sep, &ferr,
                                    &Work[0], &Work[M], &Work[2*M], NULL, &Work[3*M], &ldworkx,  &infox, 1, 1, 1, 1, 1);
                        } else {
                            mepack_double_gglyap("N", TRANSA, M, A, M, B, M, Q, M, Z, M, X, M, &scale, Work, ldwork,  &info);
                        }
                    } else {
                        if ( is == 27 ) {
                            Int infox  = 0;
                            FC_GLOBAL_(sg03ad,SG03AD)("Cont", "X", "Fact", TRANSX, "upper", &M,
                                    A, &M, B, &M, Q, &M, Z, &M, X, &M, &scale, &sep, &ferr,
                                    &Work[0], &Work[M], &Work[2*M], NULL, &Work[3*M], &ldworkx,  &infox, 1, 1, 1, 1, 1);
                        } else {
                            mepack_double_gglyap("F", TRANSA, M, A, M, B, M, Q, M, Z, M, X, M, &scale, Work, ldwork,  &info);
                        }
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

                ress += benchmark_check_X_double(M,M,X, M, Xorig, M);


            }

            times /= (double) nMAT;
            ctimes /= (double) nMAT;
            ress /= (double) nMAT;

            /* Print  */

            printf("%5d %3d %10.5e  %10.5e  %10.5e  %10.5e \n",
                   (int)  M, (int) MB, times, ctimes, ctimes/times, ress);
            fflush(stdout);
            free(A);
            free(B);
            free(Aorig);
            free(Borig);
            free(Q);
            free(Z);
            free(X);
            free(Xorig);
            free(Work);
            free(RHS);

        }
    }

    benchmark_exit();
    return (ress < sqrt(eps)*100)? 0 : -1 ;
}


