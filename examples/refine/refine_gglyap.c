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

const char * solver_name_str(int is) {
     switch(is) {
        case MEPACK_FRONTEND_SOLVER_LEVEL3:
            return ("Refine with LEVEL3 solver");
            break;
        case MEPACK_FRONTEND_SOLVER_LEVEL2:
            return ("Refine with LEVEL2 solver");
            break;
        case MEPACK_FRONTEND_SOLVER_2STAGE:
            return ("Refine with LEVEL3 two stage solver");
            break;
        case MEPACK_FRONTEND_SOLVER_DAG:
            return ("Refine with DAG solver.");
            break;
        case MEPACK_FRONTEND_SOLVER_RECURSIVE:
            return ("Refine with RECURSIVE solver.");
            break;
    }
     return "unknown";
}


void usage(char *prgmname) {
    int is;
    printf("Solve a generalized Lyapunov equation with iterative refinement.\n");
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
    printf("--alignoff, -a           Turn off the blocksize alignment.\n");
    printf("--output=path, -o path   Write the final result table to csv-like file.\n");


    printf("\nPossible Solvers\n");

    for (is = 0; is <= 3; is++) {
        printf("%2d : ", is); solver_name(is);
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
    double *convlog;
    Int M = 1024;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int MB;
    Int MAXIT = 30;
    double TAU = 0.001;
    Int RUNS = 5;
    Int is = 0;
    Int nMAT = 1;
    Int align_on = 1;
    Int ldwork;
    double alpha, beta;
    char TRANSA[20]="N";
    double res =0.0;

    int info, run;
    double te, ts;
    char *output_file = NULL;
    double times,ts2, te2;
    double ctimes;
    double ferr = 1.0;
    double eps;
    size_t mem = 1;

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
            {"solver",  required_argument, 0, 's'},
            {"alignoff", no_argument, 0, 'a'},
            {"output",  required_argument, 0, 'o'},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
            no_argument: " "
            required_argument: ":"
            optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:M:A:t:r:s:S:o:",
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
            case 'o':
                output_file = strdup ( optarg );
                break;


            default:
                /* Not sure how to get here... */
                return EXIT_FAILURE;
        }
    }

    if ( is < 1 || is > 5 ) {
        fprintf(stderr, "Solver not known. \n");
        exit(-1);
    }

    benchmark_init();
    mepack_init();

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
    csc_table_comment_printf(tab, "Solution of the Generalized Lyapunov Sylvester Equation");
    csc_table_comment_printf(tab, "Solver: %s", solver_name_str(is));

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




    mepack_tglyap_isolver_set(1);

    eps = mepack_double_epsilon();

    for (M = M_MIN ;  M <= M_MAX ; M = M + M_STEP ) {
        for (MB = MB_MIN; MB <= MB_MAX ; MB += MB_STEP) {
            /* Reset the Seed  */
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            /* Set the Block size   */
            mepack_double_tglyap_blocksize_set(MB);
            mepack_tglyap_isolver_set(1);
            mepack_tglyap_frontend_solver_set(is);
            /* Prepare  */
            times = 0;
            ctimes = 0;
            ferr = 0.0;
            res = 0.0;

            A = (double *) malloc(sizeof(double) * (M*M));
            B = (double *) malloc(sizeof(double) * (M*M));
            Aorig = (double *) malloc(sizeof(double) * (M*M));
            Borig = (double *) malloc(sizeof(double) * (M*M));
            Q = (double *) malloc(sizeof(double) * (M*M));
            Z = (double *) malloc(sizeof(double) * (M*M));



            X = (double *) malloc(sizeof(double) * (M*M));
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; B[i] = 0.0;  }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { X[i] = 0.0; }

            Xorig = (double *) malloc(sizeof(double) * (M*M));
            RHS = (double *) malloc(sizeof(double) * (M*M));

            convlog = (double *) malloc(sizeof(double) * MAXIT);
            mem = mepack_memory_frontend("DLA_GGLYAP_REFINE", "N", "N", M, M);

            Work = (double *) malloc(sizeof(double) * (mem));

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &M, &alpha, &beta, Xorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                benchmark_random_gevp_double(M, iseed, A, B, Q, Z, Aorig, Borig, align_on* MB);
                benchmark_rhs_glyap_double(TRANSA, M, Aorig, M, Borig, M, Xorig, M, RHS, M );


                te = 0.0;
                te2 = 0.0;
                for (run = -1; run < RUNS; run++) {
                    double ltau = TAU;
                    int lmaxit = MAXIT;
                    ldwork = mem;
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, RHS, &M, X, &M, 1);

                    ts = get_wtime();
                    ts2 = get_ctime();
                    mepack_double_gglyap_refine(TRANSA, "N",  M, Aorig, M, Borig, M, X, M, RHS, M,
                            A, M, B, M, Q, M, Z, M, &lmaxit, &ltau, convlog, Work, ldwork, &info);

                    if (run >= 0 ) {
                        te2 += (get_ctime()-ts2);
                        te += (get_wtime()-ts);
                    }
                    if ( info != 0) {
                        printf("Failed with INFO = %d\n", (int) info);
                    }
                }
                te = te / RUNS;
                te2 = te2/RUNS;
                times += te;
                ctimes += te2;

                res += mepack_double_residual_glyap(TRANSA, M, Aorig, M, Borig, M, X, M, RHS, M, 1.0);
                ferr += benchmark_check_X_double(M,M,X, M, Xorig, M);


            }
            times /= (double) nMAT;
            ctimes /= (double) nMAT;
            ferr /= (double) nMAT;
            res /= (double) nMAT;


            /* Print  */
            csc_table_new_row(tab);
            csc_table_set_entry_integer(tab, col_m, M);
            csc_table_set_entry_integer(tab, col_mb, MB);
            csc_table_set_entry_float(tab, col_walltime, times);
            csc_table_set_entry_float(tab, col_cputime, ctimes);
            csc_table_set_entry_float(tab, col_ratio, ctimes/times);
            csc_table_set_entry_float(tab, col_ferr, ferr);
            csc_table_set_entry_float(tab, col_res, res);
            csc_table_print_current_row(tab);


            free(convlog);
            free(A);
            free(B);
            free(X);
            free(Xorig);
            free(Work);
            free(RHS);
            free(Aorig);
            free(Borig);
            free(Q);
            free(Z);
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


