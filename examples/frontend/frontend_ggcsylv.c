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
            printf("Not Implemented\n");
            break;
        case 26:
            printf("Not Implemented\n");
            break;
    }

}

void usage(char *prgmname) {
    int is;
    printf("Solve a generalized coupled Sylvester with different solvers.\n");
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
    printf("--solver=S, -s S         Select the solver \n");
    printf("--alignoff, -a           Turn off the blocksize alignment.\n");
    printf("--changerole, -c         B and D changes places.\n");
    printf("--reuse, -R              Benchmark the factorization reuse process.\n");
    printf("--hess, -H               Solve for generalized upper Hessenberg form.\n");
    printf("--sign1=s1               Select the first sign \n");
    printf("--sign2=s2               Select the second sign \n");
    printf("--output=path, -o path   Write the final result table to csv-like file.\n");


    printf("\nPossible Solvers\n");

    for (is = 0; is <= 24; is++) {
        printf("%2d : ", is); solver_name(is);
    }
    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i, mat ;
    double *A, *B, *C, *D;
    double *Aorig, *Borig, *Corig, *Dorig;
    double *QA, *ZA, *QB, *ZB;
    double *X, *Xorig,*Work, *RHS;
    double *Y, *Yorig, *RHS2;
    Int M = 1024, N=1024, N2;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int N_MIN=1024,N_MAX=1024,N_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int NB_MIN=64,NB_MAX=64,NB_STEP=32;
    Int MB, NB;
    char *output_file = NULL;


    Int RUNS = 5;
    Int is = 0;
    Int mn_same = 0;
    Int mbnb_same = 0;
    Int nMAT = 1;
    Int align_on = 1;
    double alpha, beta;
    double sign1 = 1;
    double sign2 = 1;

    char TRANSA[20]="N";
    char TRANSB[20]="N";
    char reusestr[10];

    double scale = 1.0;
    int info, run;
    double te, ts = 0;
    size_t mem;

    double res = 0.0;
    double times,ts2 = 0, te2;
    double ctimes;
    double forward_error = 1.0;
    double eps;
    int choice;
    int changerole = 0;
    int reuse = 0, hess = 0;

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
            {"sign1",  required_argument, 0, 'S'},
            {"sign2",  required_argument, 0, 'D'},
            {"alignoff", no_argument, 0, 'a'},
            {"reuse", no_argument, 0, 'R'},
            {"hess", no_argument, 0, 'H'},
            {"output", required_argument, 0, 'o'},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
no_argument: " "
required_argument: ":"
optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:n:M:N:A:B:t:r:s:S:D:caR:Ho:",
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
                sign1 = atof(optarg);
                break;
            case 'D':
                sign2 = atof(optarg);
                break;
            case 'a':
                align_on = 0;
                break;
            case 'c':
                changerole = 1;
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
        fprintf(stderr, "Both -R and -H were selected. \n");
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
    if ( mn_same )
        csc_table_comment_printf(tab, "Cols: same as rows. %s");
    else
        csc_table_comment_printf(tab, "Cols: %d (%d:%d:%d)", (int) N , (int) N_MIN, (int) N_STEP, (int) N_MAX);
    csc_table_comment_printf(tab, "TRANSA: %s", TRANSA);
    csc_table_comment_printf(tab, "TRANSB: %s", TRANSB);
    csc_table_comment_printf(tab, "SIGN1:   %lg", sign1);
    csc_table_comment_printf(tab, "SIGN2:   %lg", sign2);
    csc_table_comment_printf(tab, "Block Alignment: %s", (align_on == 1)?"YES":"NO");
    csc_table_comment_printf(tab, "Matrix in Hessenberg Form: %d", hess);
    if ( changerole ) csc_table_comment_printf(tab, "B and D change roles");
    csc_table_comment_printf(tab, "Solution of the Generalized Coupled Sylvester Equation");

    int col_m = csc_table_add_column(tab, "M", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_n = csc_table_add_column(tab, "N", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_mb = csc_table_add_column(tab, "MB", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_nb = csc_table_add_column(tab, "NB", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_walltime = csc_table_add_column(tab, "Wall-Time", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_cputime = csc_table_add_column(tab, "CPU-Time", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_ratio = csc_table_add_column(tab, "Ratio", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_ferr  = csc_table_add_column(tab, "Forward-Error", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_res   = csc_table_add_column(tab, "Residual", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);

    csc_table_column_minwidth(tab, col_m, 5);
    csc_table_column_minwidth(tab, col_n, 5);
    csc_table_column_minwidth(tab, col_mb, 5);
    csc_table_column_minwidth(tab, col_nb, 5);
    csc_table_column_minwidth(tab, col_walltime, 12);
    csc_table_column_minwidth(tab, col_cputime, 12);
    csc_table_column_minwidth(tab, col_ratio, 12);
    csc_table_column_minwidth(tab, col_ferr, 12);
    csc_table_column_minwidth(tab, col_res, 12);
    csc_table_print_current_row(tab);

    mepack_tgcsylv_isolver_set(1);
    mepack_double_tgcsylv_blocksize_2stage_set(256);
    eps = mepack_double_epsilon();

    for (M = M_MIN, N = N_MIN; M <= M_MAX && N<=N_MAX; M = (N>=N_MAX) ? (M+M_STEP): M,  N=(N<N_MAX)?(N+N_STEP):(N_MIN)) {
        if (mn_same) {
            N = M;
        }
        for (MB = MB_MIN, NB = NB_MIN; MB <= MB_MAX && NB<=NB_MAX; MB = (NB>=NB_MAX) ? (MB+MB_STEP): MB,  NB=(NB<NB_MAX)?(NB+NB_STEP):(NB_MIN)) {
            if ( mbnb_same ) {
                NB = MB;
            }
            /* Reset the Seed  */
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            /* Set the Block size   */
            mepack_double_tgcsylv_blocksize_mb_set(MB);
            mepack_double_tgcsylv_blocksize_nb_set(NB);

            /* Prepare  */
            times = 0;
            ctimes = 0;
            forward_error = 0.0;
            res = 0.0;

            A = (double *) malloc(sizeof(double) * (M*M));
            C = (double *) malloc(sizeof(double) * (M*M));
            Aorig = (double *) malloc(sizeof(double) * (M*M));
            Corig = (double *) malloc(sizeof(double) * (M*M));
            QA = (double *) malloc(sizeof(double) * (M*M));
            ZA = (double *) malloc(sizeof(double) * (M*M));
            B = (double *) malloc(sizeof(double) * (N*N));
            D = (double *) malloc(sizeof(double) * (N*N));
            Borig  = (double *) malloc(sizeof(double) * (N*N));
            Dorig = (double *) malloc(sizeof(double) * (N*N));
            QB = (double *) malloc(sizeof(double) * (N*N));
            ZB = (double *) malloc(sizeof(double) * (N*N));

            X = (double *) malloc(sizeof(double) * (M*N));
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; C[i] = 0.0; Aorig[i] = 0.0; Corig[i]=0.0; QA[i] =0.0; ZA[i]=0.0;  }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < N*N; i+=512) { B[i] = 0.0; D[i] = 0.0;Borig[i]=0; Dorig[i]=0.0; QB[i]=0; ZB[i]=0.0;  }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*N; i+=512) { X[i] = 0.0; }

            Xorig = (double *) malloc(sizeof(double) * (M*N));
            Y = (double *) malloc(sizeof(double) * (M*N));
            Yorig = (double *) malloc(sizeof(double) * (M*N));
            RHS = (double *) malloc(sizeof(double) * (M*N));
            RHS2 = (double *) malloc(sizeof(double) * (M*N));


            if ( is < 7 )  {
                mepack_tgcsylv_isolver_set(is);
                mepack_tgcsylv_frontend_solver_set(MEPACK_FRONTEND_SOLVER_LEVEL3);
            } else if ( is < 13 ) {
                mepack_tgcsylv_isolver_set(is-6);
                mepack_tgcsylv_frontend_solver_set(MEPACK_FRONTEND_SOLVER_DAG);
            } else if ( is < 19) {
                mepack_tgcsylv_isolver_set(is-12);
                mepack_tgcsylv_frontend_solver_set(MEPACK_FRONTEND_SOLVER_2STAGE);
            } else if ( is < 24 ) {
                mepack_tgcsylv_isolver_set(is - 18);
                mepack_tgcsylv_frontend_solver_set(MEPACK_FRONTEND_SOLVER_LEVEL3);
            } else if ( is == 24 ) {
                mepack_tgcsylv_isolver_set(1);
                mepack_tgcsylv_frontend_solver_set(MEPACK_FRONTEND_SOLVER_RECURSIVE);
            }

            if ( hess ) {
                strcpy(reusestr, "Hess");
            } else if ( reuse ) {
                strcpy(reusestr, "Nofact");
            } else {
                strcpy(reusestr, "Fact");
            }
            mem = mepack_memory_frontend("DLA_GGCSYLV", reusestr, reusestr, M, N);

            Work = (double *) malloc(sizeof(double) * (mem));

            alpha = 1; beta = 1;
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &N, &alpha, &beta, Xorig, &M, 1);
            FC_GLOBAL_(dlaset,DLASET)("All", &M, &N, &alpha, &beta, Yorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                Int IDIST = 2;
                Int ldwork = mem;
                N2 = M * M;
                FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, A);
                FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, C);
                if(hess){
                    benchmark_ghess_double(M, A, QA, C, ZA, Work, ldwork);
                }
                FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, A, &M, Aorig, &M, 1);
                FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, C, &M, Corig, &M, 1);
                N2 = N * N;
                FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, B);
                FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, D);
                if(hess){
                    benchmark_ghess_double(N, B, QB, D, ZB, Work, ldwork);
                }
                FC_GLOBAL_(dlacpy,DLACPY)("All", &N, &N, B, &N, Borig, &N, 1);
                FC_GLOBAL_(dlacpy,DLACPY)("All", &N, &N, D, &N, Dorig, &N, 1);

                benchmark_rhs_ggcsylv_double(TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N,Xorig, M , Yorig, M, RHS,M , RHS2, M );


                /*-----------------------------------------------------------------------------
                 *  LEVEL 3
                 *-----------------------------------------------------------------------------*/
                te = 0.0;
                te2 = 0.0;
                for (run = -1; run < RUNS; run++) {
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS, &M, X, &M, 1);
                    FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &N, RHS2, &M, Y, &M, 1);

                    if ( run == -1 || !reuse ) {
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, Aorig, &M, A, &M, 1);
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &M, &M, Corig, &M, C, &M, 1);
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &N, &N, Borig, &N, B, &N, 1);
                        FC_GLOBAL_(dlacpy,DLACPY)("All", &N, &N, Dorig, &N, D, &N, 1);
                    }

                    ts = get_wtime();
                    ts2 = get_ctime();

                    if ( run == -1 || !reuse) {
                        if(!hess){
                            mepack_double_ggcsylv("N", "N", TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N,  QA, M, ZA, M, QB, N, ZB, N, X, M, Y, M,  &scale, Work, ldwork, &info);
                        } else{
                            mepack_double_ggcsylv("H", "H", TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N,  QA, M, ZA, M, QB, N, ZB, N, X, M, Y, M,  &scale, Work, ldwork, &info);
                        }
                    } else {
                        mepack_double_ggcsylv("F", "F", TRANSA, TRANSB, sign1, sign2,  M, N, A, M, B, N, C, M, D, N,  QA, M, ZA, M, QB, N, ZB, N, X, M, Y, M,  &scale, Work, ldwork, &info);
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

                res += mepack_double_residual_csylv(TRANSA, TRANSB, sign1, sign2, M, N, Aorig, M, Borig, N, Corig, M, Dorig, N, X, M, Y, M, RHS, M, RHS2, M, scale);
                forward_error += benchmark_check_X_double(M,N,X, M, Xorig, M);
                forward_error += benchmark_check_X_double(M,N,Y, M, Yorig, M);


            }  // nmat
            times /= (double) nMAT;
            ctimes /= (double) nMAT;
            forward_error /= (double) nMAT;
            res /= (double) nMAT;

            /* Print  */
            csc_table_new_row(tab);
            csc_table_set_entry_integer(tab, col_m, M);
            csc_table_set_entry_integer(tab, col_n, N);
            csc_table_set_entry_integer(tab, col_mb, MB);
            csc_table_set_entry_integer(tab, col_nb, NB);
            csc_table_set_entry_float(tab, col_walltime, times);
            csc_table_set_entry_float(tab, col_cputime, ctimes);
            csc_table_set_entry_float(tab, col_ratio, ctimes/times);
            csc_table_set_entry_float(tab, col_ferr, forward_error);
            csc_table_set_entry_float(tab, col_res, res);
            csc_table_print_current_row(tab);


            free(A);
            free(B);
            free(C);
            free(D);
            free(Y);
            free(Yorig);
            free(Aorig);
            free(Borig);
            free(Corig);
            free(Dorig);
            free(QA);
            free(QB);
            free(ZA);
            free(ZB);
            free(X);
            free(Xorig);
            free(Work);
            free(RHS);
            free(RHS2);

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


