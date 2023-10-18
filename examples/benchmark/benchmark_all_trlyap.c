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
#include "cscutils/table.h"


void usage(char *prgmname) {
    printf("Benchmark the different level 2 inside the level 3 solver for the standard Lyapunov equation.\n");
    printf("\n");
    printf("Usage: %s <options>\n",prgmname);
    printf("\n");
    printf("--help, -h        Display this help\n");
    printf("--rows=M , -m M   Set the number of rows to M. \n");
    printf("--rows=START:STEP:STOP   Iterate the number of rows from START to STOP with step size STEP.\n");
    printf("--mb=START:STEP:STOP     Iterate the row block size form START to STOP with step size STEP.\n");
    printf("--mb=MB                  Set the row block size.  \n");
    printf("--trans=N,T              Transpose on the matrix A \n");
    printf("--nmat=NUM               Number of different equations.\n");
    printf("--runs=R                 Number of runs per equation.\n");
    printf("--mb=START:STEP:STOP     Iterate over an interval of row blocks sizes.\n");
    printf("--mb=MB                  Set the row block size.  \n");
    printf("--res                    Compute and print the residual instead of the run time. \n");
    printf("--output=path, -o path   Write the final result table to csv-like file.\n");


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
    char *output_file = NULL;

    double alpha, beta;
    char TRANSA[20]="N";
    Int print_res = 0;
#ifdef RECSY
    double MACHINE_RECSY[12];
    Int type, j;
#endif

    double scale = 1.0;
    int info, run;
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
            {"output", required_argument, 0, 'o' },
            {0,0,0,0}
        };

        int option_index = 0;


        choice = getopt_long( argc, argv, "hm:M:A:t:r:Ro:",
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
            case 'o':
                output_file = strdup(optarg);
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
    csc_table_t *tab;
    tab = csc_table_new(0);
    int col_what = csc_table_add_column(tab, "What", CSC_TABLE_STRING, CSC_TABLE_LEFT);
    int col_m = csc_table_add_column(tab, "M", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_mb = csc_table_add_column(tab, "MB", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_recsy   = csc_table_add_column(tab, "RECSY", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_lcalign = csc_table_add_column(tab, "LC-ALIGN", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_lc      = csc_table_add_column(tab, "LOCAL-COPY", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_reorder = csc_table_add_column(tab, "REORDER", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_l2opt   = csc_table_add_column(tab, "L2-Optimized", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_unopt   = csc_table_add_column(tab, "Unoptimized", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_rec     = csc_table_add_column(tab, "Recursive", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_omp     = csc_table_add_column(tab, "OpenMP-DAG", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_l32s    = csc_table_add_column(tab, "Level3-2Stage", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);

    csc_table_comment_cmd(tab, argc, argv);
    csc_table_comment_allinfo(tab);
    csc_table_comment_printf(tab, "RUNS:  %d", (int) RUNS);
    csc_table_comment_printf(tab, "Number of Matrices: %d", (int) nMAT);
    csc_table_comment_printf(tab, "Rows: %d (%d:%d:%d)", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    csc_table_comment_printf(tab, "TRANSA: %s", TRANSA);
    csc_table_comment_printf(tab, "Solution of the Standard Lyapunov Equation");


    csc_table_column_minwidth(tab, col_what, 10);
    csc_table_column_minwidth(tab, col_m, 5);
    csc_table_column_minwidth(tab, col_mb, 5);
    csc_table_column_minwidth(tab, col_recsy, 12);
    csc_table_column_minwidth(tab, col_lcalign, 12);
    csc_table_column_minwidth(tab, col_lc, 12);
    csc_table_column_minwidth(tab, col_reorder, 12);
    csc_table_column_minwidth(tab, col_l2opt, 12);
    csc_table_column_minwidth(tab, col_unopt, 12);
    csc_table_column_minwidth(tab, col_rec, 12);
    csc_table_column_minwidth(tab, col_omp, 12);
    csc_table_column_minwidth(tab, col_l32s, 12);
    csc_table_print_current_row(tab);



    mepack_trlyap_isolver_set(1);
    if ( getenv("MEPACK_LUA_CONFIG")) {
        printf("# MEPACK_LUA_CONFIG set. Using blocksize from %s\n", getenv("MEPACK_LUA_CONFIG"));
        MB_MIN = 1;
        MB_MAX = 1;
    }


    /** printf("#  M   MB       RECSY     LC_ALIGN   LOCAL_COPY      REORDER       L2-OPT        UNOPT    RECURSIVE     OpenMP-DAG    L3-2Stage\n"); */
    for (M = M_MIN ;  M <= M_MAX ; M = M + M_STEP ) {
        for (MB = MB_MIN; MB <= MB_MAX ; MB += MB_STEP) {
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            /* Set the Block size   */
            if ( getenv("MEPACK_LUA_CONFIG")) {
                mepack_double_trlyap_blocksize_set(0);
            }
            else
            {
                mepack_double_trlyap_blocksize_set(MB);
            }


            for (i = 0; i < 9; i++) {
                times[i] = 0;
                ress[i] = 0.0;
            }

            A = (double *) malloc(sizeof(double) * (M*M));
            X = (double *) malloc(sizeof(double) * (M*M));
            Xorig = (double *) malloc(sizeof(double) * (M*M));
            RHS = (double *) malloc(sizeof(double) * (M*M));
            Work = (double *) malloc(sizeof(double) * (4*M*M));

#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0;}
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { X[i] = 0.0; }

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &M, &alpha, &beta, Xorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                benchmark_random_evp_double(M, iseed, A,  NULL, NULL, MB);
                benchmark_rhs_lyap_double (TRANSA,  M, A, M, Xorig, M, RHS, M );

#ifdef RECSY
                /*-----------------------------------------------------------------------------
                 *  RECSY
                 *-----------------------------------------------------------------------------*/
                Int infox = 0;
                type = benchmark_glyap_setup_to_type(TRANSA);
                FC_GLOBAL_(recsy_machine,RECSY_MACHINE)(MACHINE_RECSY);
                te = 0;
                j = M * M ;
                for (run = 0; run < RUNS; run++) {
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, RHS, &M, X, &M, 1);
                    ts = get_wtime();
                    FC_GLOBAL_(reclyct,RECLYCT)(&type, &scale, &M, A, &M, X, &M, &infox, MACHINE_RECSY);
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < j; i++) {
                            X[i*M+j] = X[j*M+i];
                        }
                    }
                    te += (get_wtime() - ts);
                }
                te = te / RUNS;
                times[0] += te;
                ress[0] += mepack_double_residual_lyap(TRANSA, M, A, M, X, M, RHS, M, scale);
#else
                times[0] = 0;
                ress[0] = 0;
#endif


                /*-----------------------------------------------------------------------------
                 *  Bartels-Stewart Level3
                 *-----------------------------------------------------------------------------*/
                for (is = 1; is < 9; is++) {
                    mepack_trlyap_isolver_set(is);
                    te = 0.0;
                    for (run = 0; run < RUNS; run++) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, RHS, &M, X, &M, 1);
                        if ( is < 7 ) {
                            ts = get_wtime();
                            mepack_double_trlyap_level3(TRANSA,  M, A, M, X, M, &scale, Work, &info);
                            te += (get_wtime()-ts);
                        } else if ( is == 7 ) {
                            mepack_trlyap_isolver_set(1);
                            ts = get_wtime();
                            mepack_double_trlyap_dag(TRANSA,  M, A, M, X, M, &scale, Work, &info);
                            te += (get_wtime()-ts);
                        } else if ( is == 8 ) {
                            mepack_trlyap_isolver_set(1);
                            ts = get_wtime();
                            mepack_double_trlyap_level3_2stage(TRANSA,  M, A, M, X, M, &scale, Work, &info);
                            te += (get_wtime()-ts);
                        }
                    }
                    te = te / RUNS;
                    times[is] += te;

                    ress[is] += mepack_double_residual_lyap(TRANSA, M, A, M, X, M, RHS, M, scale);
                }

            }

            for (i = 0; i < 9; i++) {
                times[i] /= (double) nMAT;
                ress[i] /= (double) nMAT;

            }

            /* Print  */
            csc_table_new_row(tab);
            csc_table_set_entry_string(tab, col_what, "Wall-Time");
            csc_table_set_entry(tab, col_m, (int)M);
            csc_table_set_entry(tab, col_mb, (int)MB);
            csc_table_set_entry(tab, col_recsy, times[0]);
            csc_table_set_entry(tab, col_lcalign, times[1]);
            csc_table_set_entry(tab, col_lc, times[2]);
            csc_table_set_entry(tab, col_reorder, times[3]);
            csc_table_set_entry(tab, col_l2opt, times[4]);
            csc_table_set_entry(tab, col_unopt, times[5]);
            csc_table_set_entry(tab, col_rec, times[6]);
            csc_table_set_entry(tab, col_omp, times[7]);
            csc_table_set_entry(tab, col_l32s, times[8]);
            csc_table_print_current_row(tab);

            if ( print_res )
            {
                csc_table_new_row(tab);
                csc_table_set_entry_string(tab, col_what, "Residual");
                csc_table_set_entry(tab, col_m, (int)M);
                csc_table_set_entry(tab, col_mb, (int)MB);
                csc_table_set_entry(tab, col_recsy, ress[0]);
                csc_table_set_entry(tab, col_lcalign, ress[1]);
                csc_table_set_entry(tab, col_lc, ress[2]);
                csc_table_set_entry(tab, col_reorder, ress[3]);
                csc_table_set_entry(tab, col_l2opt, ress[4]);
                csc_table_set_entry(tab, col_unopt, ress[5]);
                csc_table_set_entry(tab, col_rec, ress[6]);
                csc_table_set_entry(tab, col_omp, ress[7]);
                csc_table_set_entry(tab, col_l32s, ress[8]);
                csc_table_print_current_row(tab);
            }


            free(A);
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
    return 0;
}

