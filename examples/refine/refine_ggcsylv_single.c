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

void usage(char *prgmname) {
    int is;
    printf("Solve a generalized coupled Sylvester equation with iterative refinement in single precision.\n");
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
    printf("--sign1=s1               Select the sign \n");
    printf("--sign2=s2               Select the sign \n");
    printf("--alignoff, -a           Turn off the blocksize alignment.\n");
    printf("--changerole, -c         Change the role of B and D.\n");
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
    float *A, *B, *C, *D;
    float *Aorig, *Borig, *Corig, *Dorig;
    float *QA, *QB, *ZA, *ZB;
    float *Work;
    float *R, *Rorig, *RHS_R;
    float *L, *Lorig, *RHS_L;
    float *convlog;

    Int M = 1024, N=1024;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;
    Int N_MIN=1024,N_MAX=1024,N_STEP=128;
    Int MB_MIN=64,MB_MAX=64,MB_STEP=32;
    Int NB_MIN=64,NB_MAX=64,NB_STEP=32;
    Int MB, NB;
    Int align_on = 1;

    Int RUNS = 5;
    Int is = 0;
    Int mn_same = 0;
    Int mbnb_same = 0;
    Int nMAT = 1;
    Int MAXIT =30;
    float TAU = 0.001;
    char *output_file = NULL;
    float alpha, beta;
    float sign1 = 1;
    float sign2 = 1;
    char TRANSA[20]="N";
    char TRANSB[20]="N";

    int info, run;
    double te, ts;
    size_t mem;
    float res =0.0;

    double  times,ts2, te2;
    double  ctimes;
    float ferr = 1.0;
    float eps;
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
            {"sign1",  required_argument, 0, 'S'},
            {"sign2",  required_argument, 0, 'D'},
            {"alignoff", no_argument, 0, 'a'},
            {"changerole", no_argument, 0, 'c'},
            {"output",  required_argument, 0, 'o'},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
no_argument: " "
required_argument: ":"
optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:n:M:N:A:B:t:r:s:S:aco:",
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
    if ( mn_same )
        csc_table_comment_printf(tab, "Cols: same as rows. %s");
    else
        csc_table_comment_printf(tab, "Cols: %d (%d:%d:%d)", (int) N , (int) N_MIN, (int) N_STEP, (int) N_MAX);
    csc_table_comment_printf(tab, "TRANSA: %s", TRANSA);
    csc_table_comment_printf(tab, "TRANSB: %s", TRANSB);
    csc_table_comment_printf(tab, "SIGN1:   %lg", sign1);
    csc_table_comment_printf(tab, "SIGN2:   %lg", sign2);
    csc_table_comment_printf(tab, "Block Alignment: %s", (align_on == 1)?"YES":"NO");
    if ( changerole ) csc_table_comment_printf(tab, "B and D change roles");
    csc_table_comment_printf(tab, "Solution of the Generalized Coupled Sylvester Equation");
    csc_table_comment_printf(tab, "Solver: %s", solver_name_str(is));

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

    eps = mepack_single_epsilon();

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
            mepack_single_tgcsylv_blocksize_mb_set(MB);
            mepack_single_tgcsylv_blocksize_nb_set(NB);
            mepack_tgcsylv_isolver_set(1);
            mepack_tgcsylv_frontend_solver_set(is);

            /* Prepare  */
            times = 0;
            ctimes = 0;
            ferr = 0.0;
            res = 0.0;

            A = (float *) malloc(sizeof(float) * (M*M));
            C = (float *) malloc(sizeof(float) * (M*M));
            B = (float *) malloc(sizeof(float) * (N*N));
            D = (float *) malloc(sizeof(float) * (N*N));
            Aorig = (float *) malloc(sizeof(float) * (M*M));
            Corig = (float *) malloc(sizeof(float) * (M*M));
            Borig = (float *) malloc(sizeof(float) * (N*N));
            Dorig = (float *) malloc(sizeof(float) * (N*N));

            QA = (float *) malloc(sizeof(float) * (M*M));
            ZA = (float *) malloc(sizeof(float) * (M*M));
            QB = (float *) malloc(sizeof(float) * (N*N));
            ZB = (float *) malloc(sizeof(float) * (N*N));

            R = (float *) malloc(sizeof(float) * (M*N));
            L = (float *) malloc(sizeof(float) * (M*N));



            convlog = (float *) malloc(sizeof(float) * MAXIT);

#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; C[i] = 0.0;  }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < N*N; i+=512) { B[i] = 0.0; D[i] = 0.0; }
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*N; i+=512) { R[i] = 0.0; L[i] = 0.0;  }

            Rorig = (float *) malloc(sizeof(float) * (M*N));
            Lorig = (float *) malloc(sizeof(float) * (M*N));
            RHS_R = (float *) malloc(sizeof(float) * (M*N));
            RHS_L = (float *) malloc(sizeof(float) * (M*N));

            mem = mepack_memory_frontend("DLA_GGCSYLV_REFINE", "N", "N", M, N);
            Work = (float *) malloc(sizeof(float) * (mem));

            alpha = 1; beta = 1;
            FC_GLOBAL_(slaset,SLASET)("All", &M, &N, &alpha, &beta, Rorig, &M, 1);
            FC_GLOBAL_(slaset,SLASET)("All", &M, &N, &alpha, &beta, Lorig, &M, 1);


            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */

                benchmark_random_gevp_float(M, iseed, A, C, QA, ZA, Aorig, Corig, align_on*gcd(MB, NB));
                benchmark_random_gevp_float(N, iseed, B, D, QB, ZB, Borig, Dorig, align_on*gcd(MB, NB));


                benchmark_rhs_ggcsylv_float(TRANSA, TRANSB, sign1, sign2,  M, N, Aorig, M, Borig, N, Corig, M, Dorig, N, Rorig, M, Lorig, M,  RHS_R, M, RHS_L, M );

                /*-----------------------------------------------------------------------------
                 *  LEVEL 3
                 *-----------------------------------------------------------------------------*/
                                        te = 0.0;
                    te2 = 0.0;
                    for (run = -1; run < RUNS; run++) {
                        float ltau = TAU;
                        int lmaxit = MAXIT;
                        Int ldwork = mem;
                        FC_GLOBAL_(slacpy,SLACPY)("All", &M, &N, RHS_R, &M, R, &M, 1);
                        FC_GLOBAL_(slacpy,SLACPY)("All", &M, &N, RHS_L, &M, L, &M, 1);

                        ts = get_wtime();
                        ts2 = get_ctime();
                        mepack_single_ggcsylv_refine(TRANSA, TRANSB, "N", sign1, sign2,  M, N, Aorig, M, Borig, N, Corig, M, Dorig, N,  R, M, L, M, RHS_R, M, RHS_L, M,
                                A, M, B, N, C, M, D, N, QA, M, ZA, M, QB, N, ZB, N, &lmaxit, &ltau, convlog, Work, ldwork, &info);

                        if (run >= 0 ) {
                            te2 += (get_ctime()-ts2);
                            te += (get_wtime()-ts);
                        }

                    }
                    te = te / RUNS;
                    te2 = te2/RUNS;
                    times += te;
                    ctimes += te2;

                    res += mepack_single_residual_csylv(TRANSA, TRANSB, sign1, sign2, M, N, Aorig, M, Borig, N, Corig, M, Dorig, N, R, M, L, M, RHS_R, M, RHS_L, M, 1.0);
                    ferr += benchmark_check_X_float(M,N,R, M, Rorig, M);
                    ferr += benchmark_check_X_float(M,N,L, M, Lorig, M);



            }
            times /= (float) nMAT;
            ctimes /= (float) nMAT;
            ferr /= (float) nMAT;
            res /= (float) nMAT;


             /* Print  */
            csc_table_new_row(tab);
            csc_table_set_entry_integer(tab, col_m, M);
            csc_table_set_entry_integer(tab, col_n, N);
            csc_table_set_entry_integer(tab, col_mb, MB);
            csc_table_set_entry_integer(tab, col_nb, NB);
            csc_table_set_entry_float(tab, col_walltime, times);
            csc_table_set_entry_float(tab, col_cputime, ctimes);
            csc_table_set_entry_float(tab, col_ratio, ctimes/times);
            csc_table_set_entry_float(tab, col_ferr, ferr);
            csc_table_set_entry_float(tab, col_res, res);
            csc_table_print_current_row(tab);



            free(convlog);
            free(A);
            free(B);
            free(C);
            free(D);
            free(Aorig);
            free(Borig);
            free(Corig);
            free(Dorig);
            free(QA);
            free(QB);
            free(ZA);
            free(ZB);
            free(R);
            free(Rorig);
            free(RHS_R);
            free(L);
            free(Lorig);
            free(RHS_L);

            free(Work);
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


