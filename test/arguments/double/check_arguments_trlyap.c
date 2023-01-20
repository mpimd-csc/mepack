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
#include <complex.h>
#include <math.h>
#include "proto.h"
#include <string.h>
typedef void (*trlyap_function) (const char *TRANSA, Int *M, double * A, Int *LDA,double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1);


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

#define TRLYAP_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9)  do { cname = trlyap_functions[i].name; global_state =0;  \
            ign_wrng = trlyap_functions[i].ign_wrng_name, count_all++;  \
            trlyap(A1,A2,A3,A4,A5,A6,A7,A8,A9,1); \
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
    Int M = 0, LDA=0, LDX=0, INFO=0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    trlyap_function trlyap;
    struct function_name {
        trlyap_function fn;
        char * name;
        int chk_l3;
        int chk_l3_2s;
        int ign_wrng_name;
    };

    struct function_name trlyap_functions[] = {
        { &(FC_GLOBAL_(dla_trlyap_l3,DLA_TRLYAP_L3)), "DLA_TRLYAP_L3", 1, 0,0 },
        { &(FC_GLOBAL_(dla_trlyap_l2,DLA_TRLYAP_L2)), "DLA_TRLYAP_L2", 0, 0 ,0 },
        { &(FC_GLOBAL_(dla_trlyap_l2_opt,DLA_TRLYAP_L2_OPT)), "DLA_TRLYAP_L2_OPT", 0,0, 0 },
        { &(FC_GLOBAL_(dla_trlyap_dag,DLA_TRLYAP_DAG)), "DLA_TRLYAP_DAG", 1, 0, 0 },
        { &(FC_GLOBAL_(dla_trlyap_l3_2s,DLA_TRLYAP_L3_2S)), "DLA_TRLYAP_L3_2S",1 , 1 ,0},
        { &(FC_GLOBAL_(dla_trlyap_recursive,DLA_TRLYAP_RECURSIVE)), "DLA_TRLYAP_RECURSIVE", 0, 0 , 0},
        { NULL, NULL, 0,0, 0}
    };


    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( trlyap_functions[i].fn != NULL ) {
        printf("Check %s\n", trlyap_functions[i].name);
        trlyap = trlyap_functions[i].fn;

        /* Check detection of argument 1   */
        TRLYAP_CALL(1, "X", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
        /* Check detection of argument 2   */
        M = -1;  TRLYAP_CALL(2, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0;
        /* Check detection of argument 4   */
        LDA = 0;  TRLYAP_CALL(4, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
        M=2; LDA = 1;  TRLYAP_CALL(4, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
        LDA = 5; M = 0;

        /* Check detection of argument 6   */
        LDX = 0;  TRLYAP_CALL(6, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
        M=2; LDX = 1;  TRLYAP_CALL(6, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
        LDX = 1; M = 0;


        if ( trlyap_functions[i].chk_l3 ) {
            mepack_double_trlyap_blocksize_set(1);
            TRLYAP_CALL(30, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_trlyap_blocksize_set(2);
        }
        if ( trlyap_functions[i].chk_l3_2s ) {
            mepack_double_trlyap_blocksize_set(4);
            mepack_double_trlyap_blocksize_2stage_set(2);
            TRLYAP_CALL(33, "N", &M, NULL, &LDA, NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_trlyap_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

