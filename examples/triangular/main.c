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

#define MAX_SOLVER 40

#ifdef SINGLE_PRECISION
#define MEPACK_PREFIX(x) mepack_single_ ## x
#define FLOAT float
#else
#define MEPACK_PREFIX(x) mepack_double_ ## x
#define FLOAT double
#endif

benchmark_context_t *benchmark_create(void);
const char *equation_name(void);
int have_solver(int is);
int is_sym(void);

const char * solver_name_str(int is) {
    switch(is) {
        case 1:
            return "LEVEL3 - LEVEL2: LOCAL COPY fixed maximum size with alignment";
            break;
        case 2:
            return "LEVEL3 - LEVEL2: LOCAL COPY fixed maximum size";
            break;
        case 3:
            return "LEVEL3 - LEVEL2: REORDERED  column first solution";
            break;
        case 4:
            return "LEVEL3 - LEVEL2: BLAS LEVEL-2 call";
            break;
        case 5:
            return "LEVEL3 - LEVEL2: Unoptimized";
            break;
        case 6:
            return "LEVEL3 - LEVEL2: RECURSIVE BLOCKING";
            break;
        case 7:
            return "DAG - LEVEL2: LOCAL COPY fixed maximum size with alignment";
            break;
        case 8:
            return "DAG - LEVEL2: LOCAL COPY fixed maximum size";
            break;
        case 9:
            return "DAG - LEVEL2: REORDERED  column first solution";
            break;
        case 10:
            return "DAG - LEVEL2: BLAS LEVEL-2 call";
            break;
        case 11:
            return "DAG - LEVEL2: Unoptimized";
            break;
        case 12:
            return "DAG - LEVEL2: RECURSIVE BLOCKING";
            break;
        case 13:
            return "LEVEL3 - DAG - LEVEL2: LOCAL COPY fixed maximum size with alignment";
            break;
        case 14:
            return "LEVEL3 - DAG - LEVEL2: LOCAL COPY fixed maximum size";
            break;
        case 15:
            return "LEVEL3 - DAG - LEVEL2: REORDERED  column first solution";
            break;
        case 16:
            return "LEVEL3 - DAG - LEVEL2: BLAS LEVEL-2 call";
            break;
        case 17:
            return "LEVEL3 - DAG - LEVEL2: Unoptimized";
            break;
        case 18:
            return "LEVEL3 - DAG - LEVEL2: RECURSIVE BLOCKING";
            break;
        case 19:
            return "LEVEL2: LOCAL COPY fixed maximum size with alignment";
            break;
        case 20:
            return "LEVEL2: LOCAL COPY fixed maximum size";
            break;
        case 21:
            return "LEVEL2: REORDERED  column first solution";
            break;
        case 22:
            return "LEVEL2: BLAS LEVEL-2 call";
            break;
        case 23:
            return "LEVEL2: Unoptimized";
            break;
        case 24:
            return "RECURSIVE BLOCKING";
            break;
        case 25:
            return "RECSY";
            break;
        case 26:
            return "RECSY PARALLEL ";
            break;
        case 27:
            return "Gardiner-Laub Algorithm";
            break;
        case 28:
            return "Gardiner-Laub Algorithm (Algorithm 705 implementation)";
            break;
        case 29:
            return "LAPACK";
            break;
        case 30:
            return "SLICOT";
            break;
    }

    return "";



}

void usage(char *prgmname) {
    int is;
    printf("Solve a %s.\n", equation_name());
    printf("\n");
    printf("Usage: %s <options>\n", prgmname);
    printf("\n");
    printf("--help, -h        Display this help\n");
    printf("--rows=M , -m M   Set the number of rows to M. \n");
    printf("--mb=START:STEP:STOP       Iterate the row block size form START to STOP with step size STEP.\n");
    printf("--mb=MB                    Set the row block size.  \n");
    printf("--transa=N,T, --trans=N,T  Transpose on the matrix pair (A,C) \n");
    if (!is_sym()) {
        printf("--cols=N , -n N   Set the number of cols to N. \n");
        printf("--cols=START:STEP:STOP   Iterate the number of columns from START to STOP with step size STEP.\n");
        printf("--cols=-                 Use the number of rows as number of cols.\n");
        printf("--nb=START:STEP:STOP     Iterate the column column block size from START to STOP with step size STEP.\n");
        printf("--nb=NB                  Set the column block size\n");
        printf("--nb=-                   Set the column block size to the row block size.\n");
        printf("--transb=N,T             Transpose on the matrix pair (B,D) \n");
    }
    printf("--nmat=NUM               Number of different equations.\n");
    printf("--runs=R                 Number of runs per equation.\n");
    printf("--solver=S, -s S         Select the solver \n");
    printf("--sign1=s1               Select the sign in the first equation. (CSYLV only)\n");
    printf("--sign2=s2               Select the sign in the second equation.(CSYLV only)\n");
    printf("--sign=s                 Select the sign between the terms of the equation. (Sylvester only)\n");
    printf("--alignoff, -a           Turn off the blocksize alignment.\n");
    printf("--changerole, -c         Change the role of B and D. (Sylvester only)\n");
    printf("--output file, -o file   Output the result table.\n");
    printf("\nPossible Solvers\n");

    for (is = 0; is <= MAX_SOLVER; is++) {
        if ( have_solver(is))
            printf("%2d : %s\n", is, solver_name_str(is));
    }
    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int mat ;
    char *output_file = NULL;

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

    FLOAT sign1 = 1;
    FLOAT sign2 = 1;
    char TRANSA[20]="N";
    char TRANSB[20]="N";
    csc_table_t *tab;

    int run;
    double te, ts = 0;

    double times,ts2 = 0, te2;
    double ctimes;
    double forwarderr = 1.0;
    double relres = 0.0;
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
            {"trans",   required_argument, 0, 'A'},
            {"transb",  required_argument, 0, 'B'},
            {"nmat",    required_argument, 0, 't'},
            {"runs",    required_argument, 0, 'r'},
            {"solver",  required_argument, 0, 's'},
            {"sign1",  required_argument, 0, 'S'},
            {"sign2",  required_argument, 0, 'D'},
            {"sign",   required_argument, 0, 'S'},
            {"alignoff", no_argument, 0, 'a'},
            {"changerole", no_argument, 0, 'c'},
            {"output", required_argument, 0, 'o'},
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
                output_file = strdup(optarg);
                break;
            default:
                /* Not sure how to get here... */
                return EXIT_FAILURE;
        }
    }

    if ( ! have_solver(is) ) {
        fprintf(stderr, "Solver not known or not supported. \n");
        exit(0);
    }

    benchmark_init();
    mepack_init();


    /*-----------------------------------------------------------------------------
     *  Output Configuration
     *-----------------------------------------------------------------------------*/
    tab = csc_table_new(0);
    csc_table_comment_cmd(tab, argc, argv);
    csc_table_comment_allinfo(tab);
    csc_table_comment_printf(tab, "RUNS:  %d", (int) RUNS);
    csc_table_comment_printf(tab, "Number of Matrices: %d", (int) nMAT);
    csc_table_comment_printf(tab, "Rows: %d (%d:%d:%d)", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);
    if (!is_sym()) {
    if ( mn_same )
        csc_table_comment_printf(tab, "Cols: same as rows. %s");
    else
        csc_table_comment_printf(tab, "Cols: %d (%d:%d:%d)", (int) N , (int) N_MIN, (int) N_STEP, (int) N_MAX);
    }
    csc_table_comment_printf(tab, "TRANSA: %s", TRANSA);
    if ( !is_sym()) csc_table_comment_printf(tab, "TRANSB: %s", TRANSB);
    csc_table_comment_printf(tab, "SIGN1:   %lg", sign1);
    csc_table_comment_printf(tab, "SIGN2:   %lg", sign2);
    csc_table_comment_printf(tab, "Block Alignment: %s", (align_on == 1)?"YES":"NO");
    if ( changerole ) csc_table_comment_printf(tab, " B and D change roles.");
    csc_table_comment_printf(tab, "Solver: %s", solver_name_str(is));
    csc_table_comment_printf(tab, "");
    csc_table_comment_printf(tab, "Solution of the %s:", equation_name());

    int col_m = csc_table_add_column(tab, "M", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_n = csc_table_add_column(tab, "N", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_mb = csc_table_add_column(tab, "MB", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_nb = csc_table_add_column(tab, "NB", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT);
    int col_walltime = csc_table_add_column(tab, "Wall-Time", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_cputime  = csc_table_add_column(tab, "CPU-Time", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_ratio    = csc_table_add_column(tab, "Ratio", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_forward   = csc_table_add_column(tab, "Forward-Err.", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);
    int col_relres   = csc_table_add_column(tab, "Rel. Residual", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT);

    csc_table_column_minwidth(tab, col_m, 5);
    csc_table_column_minwidth(tab, col_n, 5);
    csc_table_column_minwidth(tab, col_mb, 5);
    csc_table_column_minwidth(tab, col_nb, 5);
    csc_table_column_minwidth(tab, col_walltime, 12);
    csc_table_column_minwidth(tab, col_cputime, 12);
    csc_table_column_minwidth(tab, col_ratio, 12);
    csc_table_column_minwidth(tab, col_forward, 12);
    csc_table_column_minwidth(tab, col_relres, 12);
    csc_table_print_current_row(tab);


    mepack_tgcsylv_isolver_set(1);
    eps = MEPACK_PREFIX(epsilon)();


    benchmark_context_t *ctx = benchmark_create();

    if ( is_sym()) {
        mn_same = 1;
        mbnb_same = 1;
        N_MIN = M_MIN;
        N_MAX = M_MAX;
        N_STEP = M_STEP;
        NB_MIN = MB_MIN;
        NB_MAX = MB_MAX;
        NB_STEP = MB_STEP;
    }


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

            /* Prepare  */
            times = 0;
            ctimes = 0;
            forwarderr = 0.0;
            relres = 0.0;

            benchmark_context_init(ctx,M, N, is, MB, NB, 256, TRANSA, TRANSB, sign1, sign2);
            benchmark_context_setup_rhs(ctx);

            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                benchmark_context_setup_problem(ctx, iseed, align_on*gcd(MB, NB), changerole);

                te = 0.0;
                te2 = 0.0;
                for (run = -1; run < RUNS; run++) {
                    benchmark_context_iter_prepare(ctx);

                    ts = get_wtime();
                    ts2 = get_ctime();
                    benchmark_context_iter_solve(ctx);
                    if (run >= 0 ) {
                        te2 += (get_ctime()-ts2);
                        te += (get_wtime()-ts);
                    }

                }
                te = te / RUNS;
                te2 = te2/RUNS;
                times += te;
                ctimes += te2;

                benchmark_context_check(ctx, &forwarderr, &relres);
            }
            times /= (double) nMAT;
            ctimes /= (double) nMAT;
            forwarderr /= (double) nMAT;
            relres /= (double) nMAT;

            /* Print  */

            csc_table_new_row(tab);
            csc_table_set_entry(tab, col_m, (int)M);
            if(!is_sym()) csc_table_set_entry(tab, col_n, (int)N);
            csc_table_set_entry(tab, col_mb, (int)MB);
            if(!is_sym()) csc_table_set_entry(tab, col_nb, (int)NB);
            csc_table_set_entry(tab, col_walltime, (double) times);
            csc_table_set_entry(tab, col_cputime, (double) ctimes);
            csc_table_set_entry(tab, col_ratio, (double) ctimes/times);
            csc_table_set_entry(tab, col_forward, forwarderr);
            csc_table_set_entry(tab, col_relres, relres);
            csc_table_print_current_row(tab);


            fflush(stdout);
            benchmark_context_clean(ctx);
        }
    }

    benchmark_context_free(ctx);

    if ( output_file) {
        FILE *fp = fopen(output_file, "w");
        csc_table_print_ascii(fp, tab, " ");
        fclose(fp);
        free(output_file);
    }
    csc_table_destroy(tab);
    benchmark_exit();
    return ( relres < eps*1000)? 0 : -1 ;

}


