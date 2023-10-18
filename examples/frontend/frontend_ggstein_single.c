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
    printf("Solve a generalized Stein equation with different solvers in single precision.\n");
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
    printf("--hess, -H               Solve for a generalized Hessenberg form (H,T).\n");
    printf("--output=path, -o path   Write the final result table to csv-like file.\n");



    printf("\nPossible Solvers\n");
    for (is = 0; is <= 24; is++) {
        printf("%2d : ", is);
        solver_name(is);
    }

    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i,mat ;
    float *A, *B;
    float *Aorig, *Borig;
    float *Q, *Z;
    float *X, *Xorig,*Work, *RHS;
    Int M = 1024, N2;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int MB;
    char *output_file = NULL;


    Int RUNS = 5;
    Int is = 0;
    Int nMAT = 1;
    Int align_on = 1;
    size_t mem = 0 ;

    float alpha, beta;
    char TRANSA[20]="N";
    char reusestr[10];
    float eps;

    float res = 0.0;
    float scale = 1.0;
    int info, run;
    double te, ts = 0;

    double times,ts2 = 0, te2;
    double ctimes;
    float forward_error = 1.0;

    int choice;
    int reuse = 0, hess = 0;

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
            {"hess", no_argument, 0, 'H'},
            {"output", required_argument, 0, 'o'},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
no_argument: " "
required_argument: ":"
optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:M:A:t:r:s:S:R:Ho:",
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
                align_on = 0;
                break;
            case 'R':
                reuse = 1;
                break;
            case 'H':
                hess = 1;
                break;
            case 'o':
                output_file = strdup(optarg);
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

    if ( reuse == 1 && hess == 1) {
        fprintf(stderr, "-R and -H can not be selected together. \n");
        exit(-1);
    }

    mepack_init();
    benchmark_init();
    /*-----------------------------------------------------------------------------
     *  Output Configuration
     *-----------------------------------------------------------------------------*/
    csc_table_t *tab = csc_table_new(0);
    csc_table_comment_cmd(tab, argc, argv);
    csc_table_comment_allinfo(tab);
    csc_table_comment_printf(tab, "RUNS:  %d", (int) RUNS);
    csc_table_comment_printf(tab, "Number of Matrices: %d", (int) nMAT);
    csc_table_comment_printf(tab, "Rows: %d (%d:%d:%d)", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    csc_table_comment_printf(tab, "TRANSA: %s", TRANSA);
    csc_table_comment_printf(tab, "Block Alignment: %s", (align_on == 1)?"YES":"NO");
    csc_table_comment_printf(tab, "Matrix in Hessenberg Form: %d", hess);
    csc_table_comment_printf(tab, "Solution of the Generalized Stein Equation");

    int col_m = csc_table_add_column(tab, "M", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_mb = csc_table_add_column(tab, "MB", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_walltime = csc_table_add_column(tab, "Wall-Time", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_cputime = csc_table_add_column(tab, "CPU-Time", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_ratio = csc_table_add_column(tab, "Ratio", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_ferr  = csc_table_add_column(tab, "Forward-Error", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_res   = csc_table_add_column(tab, "Residual", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);

    csc_table_column_minwidth(tab, col_m, 5);
    csc_table_column_minwidth(tab, col_mb, 5);
    csc_table_column_minwidth(tab, col_walltime, 12);
    csc_table_column_minwidth(tab, col_cputime, 12);
    csc_table_column_minwidth(tab, col_ratio, 12);
    csc_table_column_minwidth(tab, col_ferr, 12);
    csc_table_column_minwidth(tab, col_res, 12);
    csc_table_print_current_row(tab);


    mepack_tgstein_isolver_set(1);
    mepack_single_tgstein_blocksize_2stage_set(256);
    eps = mepack_single_epsilon();

    for (M = M_MIN ;  M <= M_MAX ; M = M + M_STEP ) {
        for (MB = MB_MIN; MB <= MB_MAX ; MB += MB_STEP) {
            /* Reset the Seed  */
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;


            /* Set the Block size   */
            mepack_single_tgstein_blocksize_set(MB);

            /* Prepare  */
            times = 0;
            ctimes = 0;
            forward_error = 0.0;
            res = 0.0;

            A = (float *) malloc(sizeof(float) * (M*M));
            B = (float *) malloc(sizeof(float) * (M*M));
            Q = (float *) malloc(sizeof(float) * (M*M));
            Z = (float *) malloc(sizeof(float) * (M*M));
            Aorig = (float *) malloc(sizeof(float) * (M*M));
            Borig = (float *) malloc(sizeof(float) * (M*M));
            X = (float *) malloc(sizeof(float) * (M*M));
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; B[i] = 0.0; Q[i] =0.0; Z[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { X[i] = 0.0; }

            Xorig = (float *) malloc(sizeof(float) * (M*M));
            RHS = (float *) malloc(sizeof(float) * (M*M));


            if ( is < 7 )  {
                mepack_tgstein_isolver_set(is);
                mepack_tgstein_frontend_solver_set(MEPACK_FRONTEND_SOLVER_LEVEL3);
            } else if ( is > 6 && is < 13 )  {
                mepack_tgstein_isolver_set(is-6);
                mepack_tgstein_frontend_solver_set(MEPACK_FRONTEND_SOLVER_DAG);
            } else if ( is > 12 && is < 19 )  {
                mepack_tgstein_isolver_set(is-12);
                mepack_tgstein_frontend_solver_set(MEPACK_FRONTEND_SOLVER_2STAGE);
            } else if ( is == 24  ) {
                mepack_tgstein_isolver_set(1);
                mepack_tgstein_frontend_solver_set(MEPACK_FRONTEND_SOLVER_RECURSIVE);
            }


            if ( hess ) {
                strcpy(reusestr, "Hess");
            } else if ( reuse ) {
                strcpy(reusestr, "Nofact");
            } else {
                strcpy(reusestr, "Fact");
            }
            mem = mepack_memory_frontend("SLA_GGSTEIN", reusestr, reusestr, M, M);
            Work = (float *) malloc(sizeof(float) * (mem));

            alpha = 1; beta = 1;
            FC_GLOBAL_(slaset,SLASET)("All", &M, &M, &alpha, &beta, Xorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                // benchmark_random_gevp_float(M, iseed, A, B, NULL, NULL, NULL, NULL, align_on*MB);
                Int IDIST = 2;
                Int ldwork = mem;
                N2 = M * M;
                FC_GLOBAL(slarnv,SLARNV)(&IDIST, iseed, &N2, A);
                FC_GLOBAL(slarnv,SLARNV)(&IDIST, iseed, &N2, B);
                if(hess){
                    benchmark_ghess_single(M, A, Q, B, Z, Work, ldwork);
                }
                FC_GLOBAL_(slacpy,SLACPY)("All", &M, &M, A, &M, Aorig, &M, 1);
                FC_GLOBAL_(slacpy,SLACPY)("All", &M, &M, B, &M, Borig, &M, 1);
                benchmark_rhs_gstein_float(TRANSA, M, A, M, B, M, Xorig, M, RHS, M );

                te = 0.0;
                te2 = 0.0;
                for (run = -1; run < RUNS; run++) {
                    FC_GLOBAL_(slacpy,SLACPY)("All", &M, &M, RHS, &M, X, &M, 1);
                    if ( run == -1 || !reuse ) {
                        FC_GLOBAL_(slacpy,SLACPY)("All", &M, &M, Aorig, &M, A, &M, 1);
                        FC_GLOBAL_(slacpy,SLACPY)("All", &M, &M, Borig, &M, B, &M, 1);
                    }

                    ts = get_wtime();
                    ts2 = get_ctime();

                    if ( run == -1 || !reuse) {
                        if ( hess ) {
                            mepack_single_ggstein("H", TRANSA, M, A, M, B, M, Q, M, Z, M, X, M, &scale, Work, ldwork,  &info);
                        } else {
                            mepack_single_ggstein("N", TRANSA, M, A, M, B, M, Q, M, Z, M, X, M, &scale, Work, ldwork,  &info);
                        }

                    } else {
                        mepack_single_ggstein("F", TRANSA, M, A, M, B, M, Q, M, Z, M, X, M, &scale, Work, ldwork,  &info);
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

                res += mepack_single_residual_gstein(TRANSA, M, Aorig, M, Borig, M, X, M, RHS, M, scale);
                forward_error += benchmark_check_X_float(M,M,X, M, Xorig, M);

            }

            times /= (float) nMAT;
            ctimes /= (float) nMAT;
            forward_error /= (float) nMAT;
            res /= (float) nMAT;

            /* Print  */
            csc_table_new_row(tab);
            csc_table_set_entry_integer(tab, col_m, M);
            csc_table_set_entry_integer(tab, col_mb, MB);
            csc_table_set_entry_float(tab, col_walltime, times);
            csc_table_set_entry_float(tab, col_cputime, ctimes);
            csc_table_set_entry_float(tab, col_ratio, ctimes/times);
            csc_table_set_entry_float(tab, col_ferr, forward_error);
            csc_table_set_entry_float(tab, col_res, res);
            csc_table_print_current_row(tab);

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
    if ( output_file) {
        FILE *fp = fopen(output_file, "w");
        csc_table_print_ascii(fp, tab, " ");
        fclose(fp);
        free(output_file);
    }
    csc_table_destroy(tab);


    benchmark_exit();
    return ( res < eps*100)? 0 : -1 ;
}


