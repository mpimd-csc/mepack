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

void usage(void) {
    printf("Generate random matrix pairs and their generalized Schur decompositions for the usage with the triangular solver examples.\n\n") ;

    printf("Usage: <prgm> <options>\n");
    printf("\n");
    printf("--help, -h        Display this help\n");
    printf("--rows=M , -m M   Set the number of rows to M. \n");
    printf("--rows=START:STEP:STOP   Iterate the number of rows from START to STOP with step size STEP.\n");
    printf("--nmat=NUM               Number of different equations.\n");
    return;
}


int main(int argc, char **argv)
{
    Int iseed[4]={1,1,1,9};
    Int i,mat ;
    double *A, *B;
    Int M = 1024;
    Int M_MIN=1024,M_MAX=1024,M_STEP=128;

    Int nMAT = 1;

    int choice;


    while (1)
    {
        static struct option long_options[] =
        {
            /* Use flags like so:
            {"verbose",    no_argument,    &verbose_flag, 'V'}*/
            /* Argument styles: no_argument, required_argument, optional_argument */
            {"help",    no_argument,    0,    'h'},
            {"rows",    required_argument, 0, 'm'},
            {"nmat",    required_argument, 0, 't'},
            {0,0,0,0}
        };

        int option_index = 0;

        /* Argument parameters:
            no_argument: " "
            required_argument: ":"
            optional_argument: "::" */

        choice = getopt_long( argc, argv, "hm:t:",
                    long_options, &option_index);

        if (choice == -1)
            break;

        switch( choice )
        {
            case 'h':
                usage();
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
            case 't':
                nMAT = atoi(optarg);
                break;
            default:
                /* Not sure how to get here... */
                return EXIT_FAILURE;
        }
    }

    benchmark_init();
    /*-----------------------------------------------------------------------------
     *  Output Configuration
     *-----------------------------------------------------------------------------*/
    printf("# Command Line: ");
    for (i = 1; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("# Number of Matrices: %d\n", (int) nMAT);
    printf("# Rows: %d (%d:%d:%d)\n", (int) M , (int) M_MIN, (int) M_STEP, (int) M_MAX);


    for (M = M_MIN ;  M <= M_MAX ; M = M + M_STEP ) {
            /* Reset the Seed  */
            iseed[0] = 1;
            iseed[1] = 1;
            iseed[2] = 1;
            iseed[3] = 1;

            /* Prepare  */

            A = (double *) malloc(sizeof(double) * (M*M));
            B = (double *) malloc(sizeof(double) * (M*M));
#pragma omp parallel for schedule(static,1)
            for (i = 0; i < M*M; i+=512) { A[i] = 0.0; B[i] = 0.0;  }

            for (mat = 0; mat < nMAT; mat++) {
                /* Setup the Problem  */
                benchmark_random_gevp_double(M, iseed, A, B, NULL, NULL, NULL, NULL, 0);

            }
            fflush(stdout);
            free(A);
            free(B);
    }

    benchmark_exit();
    return 0;
}


