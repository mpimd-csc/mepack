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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "proto.h"

typedef void (*tgstein_function) (const char *TRANSA, Int *M, double * A, Int *LDA, double *, Int *LDB, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);


Int global_state = 0;
char *cname = NULL;
int ign_wrng;

void err_handler(const char * str, Int Info)
{
    if ( strcmp(str,cname) == 0 || ign_wrng )
        global_state = Info;
    else {
        printf("-- Wrong Name in error handler.\n");
        global_state = -Info;
    }
}

#define TGSTEIN_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11)  do { cname = tgstein_functions[i].name; global_state =0;  \
            ign_wrng = tgstein_functions[i].ign_wrng_name, count_all++;  \
            tgstein(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,1); \
            if ( global_state == D ) {\
                ok &= 1; count_ok ++;  \
                printf("-- ARGUMENT %d OK.\n",D); \
            } else { \
                printf("-- ARGUMENT %d FAILED (info = %d).\n",D, (int) global_state); \
                ok &=0; count_failed++;   \
            } \
        } while(0);

int main(int argc, char **argv)
{
    double WORK = 0 , SCALE = 0;
    Int M = 0, LDA=0, LDB=0, LDX=0, INFO=0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    tgstein_function tgstein;
    struct function_name {
        tgstein_function fn;
        char * name;
        int chk_l3;
        int chk_l3_2s;
        int ign_wrng_name;
    };

    struct function_name tgstein_functions[] = {
        { &(FC_GLOBAL_(dla_tgstein_l3,DLA_TGSTEIN_L3)), "DLA_TGSTEIN_L3", 1, 0,0 },
        { &(FC_GLOBAL_(dla_tgstein_l2,DLA_TGSTEIN_L2)), "DLA_TGSTEIN_L2", 0, 0 ,0 },
        { &(FC_GLOBAL_(dla_tgstein_dag,DLA_TGSTEIN_DAG)), "DLA_TGSTEIN_DAG", 1, 0, 0 },
        { &(FC_GLOBAL_(dla_tgstein_l3_2s,DLA_TGSTEIN_L3_2S)), "DLA_TGSTEIN_L3_2S",1 , 1 ,0},
        { &(FC_GLOBAL_(dla_tgstein_recursive,DLA_TGSTEIN_RECURSIVE)), "DLA_TGSTEIN_RECURSIVE", 0, 0 , 0},
        { NULL, NULL, 0,0, 0}
    };


    argc = argc;
    (void) argv;


    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( tgstein_functions[i].fn != NULL ) {
        printf("Check %s\n", tgstein_functions[i].name);
        tgstein = tgstein_functions[i].fn;

        /* Check detection of argument 1   */
        TGSTEIN_CALL(1, "X", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        /* Check detection of argument 2   */
        M = -1;  TGSTEIN_CALL(2, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0;
        /* Check detection of argument 4   */
        LDA = 0;  TGSTEIN_CALL(4, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M=2; LDA = 1;  TGSTEIN_CALL(4, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
        LDA = 5; M = 0;

        /* Check detection of argument 6   */
        LDB = 0;  TGSTEIN_CALL(6, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M=2; LDB = 1;  TGSTEIN_CALL(6, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
        LDB = 5; M = 0;

        /* Check detection of argument 8   */
        LDX = 0;  TGSTEIN_CALL(8, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M=2; LDX = 1;  TGSTEIN_CALL(8, "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        LDX = 1; M = 0;


        if ( tgstein_functions[i].chk_l3 ) {
            mepack_double_tgstein_blocksize_set(1);
            TGSTEIN_CALL(30, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_tgstein_blocksize_set(2);
        }
        if ( tgstein_functions[i].chk_l3_2s ) {
            mepack_double_tgstein_blocksize_set(4);
            mepack_double_tgstein_blocksize_2stage_set(2);
            TGSTEIN_CALL(33, "N", &M, NULL, &LDA,NULL, &LDB,  NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_tgstein_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

