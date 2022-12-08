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
#include <complex.h>
#include <math.h>
#include "proto.h"
#include <string.h>
typedef void (*trsylv_function) (const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


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

#define TRSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)  do { cname = trsylv_functions[i].name; global_state =0 ;  \
            ign_wrng = trsylv_functions[i].ign_wrng_name, count_all++;  \
            trsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,1,1); \
            if ( global_state == D ) {\
                ok &= 1; count_ok ++;  \
                printf("-- ARGUMENT %d OK.\n",D); \
            } else { \
                printf("-- ARGUMENT %d FAILED.\n",D); \
                ok &=0; count_failed++;   \
            } \
        } while(0);

int main(int argc, char **argv)
{
    double WORK = 0 , sign = 0, SCALE = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, LDX=0, INFO=0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    trsylv_function trsylv;
    struct function_name {
        trsylv_function fn;
        char * name;
        int chk_l3;
        int chk_l3_2s;
        int ign_wrng_name;
    };

    struct function_name trsylv_functions[] = {
        { &(FC_GLOBAL_(dla_trsylv_l3,DLA_TRSYLV_L3)), "DLA_TRSYLV_L3", 1, 0,0 },
        { &(FC_GLOBAL_(dla_trsylv_l3_unopt,DLA_TRSYLV_L3_UNOPT)), "DLA_TRSYLV_L3_UNOPT", 1, 0,0  },
        { &(FC_GLOBAL_(dla_trsylv_l2,DLA_TRSYLV_L2)), "DLA_TRSYLV_L2", 0, 0 ,0 },
        { &(FC_GLOBAL_(dla_trsylv_l2_reorder,DLA_TRSYLV_L2_REORDER)), "DLA_TRSYLV_L2_REORDER", 0, 0,0  },
        { &(FC_GLOBAL_(dla_trsylv_l2_local_copy,DLA_TRSYLV_L2_LOCAL_COPY)), "DLA_TRSYLV_L2_LOCAL_COPY" ,0 ,0,0  },
        { &(FC_GLOBAL_(dla_trsylv_l2_local_copy_32,DLA_TRSYLV_L2_LOCAL_COPY_32)), "DLA_TRSYLV_L2_LOCAL_COPY_32", 0,0 , 1},
        { &(FC_GLOBAL_(dla_trsylv_l2_local_copy_64,DLA_TRSYLV_L2_LOCAL_COPY_64)), "DLA_TRSYLV_L2_LOCAL_COPY_64", 0,0 , 1},
        { &(FC_GLOBAL_(dla_trsylv_l2_local_copy_96,DLA_TRSYLV_L2_LOCAL_COPY_96)), "DLA_TRSYLV_L2_LOCAL_COPY_96", 0,0 , 1},
        { &(FC_GLOBAL_(dla_trsylv_l2_local_copy_128,DLA_TRSYLV_L2_LOCAL_COPY_128)), "DLA_TRSYLV_L2_LOCAL_COPY_128", 0,0,1 },
        { &(FC_GLOBAL_(dla_trsylv_l2_unopt,DLA_TRSYLV_L2_UNOPT)), "DLA_TRSYLV_L2_UNOPT", 0,0, 0 },
        { &(FC_GLOBAL_(dla_trsylv_dag,DLA_TRSYLV_DAG)), "DLA_TRSYLV_DAG", 1, 0, 0 },
        { &(FC_GLOBAL_(dla_trsylv_l3_2s,DLA_TRSYLV_L3_2S)), "DLA_TRSYLV_L3_2S",1 , 1 ,0},
        { &(FC_GLOBAL_(dla_trsylv_recursive,DLA_TRSYLV_RECURSIVE)), "DLA_TRSYLV_RECURSIVE", 0, 0 , 0},
        { NULL, NULL, 0,0, 0}
    };


    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( trsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", trsylv_functions[i].name);
        trsylv = trsylv_functions[i].fn;

        /* Check detection of argument 1   */
        TRSYLV_CALL(1, "X", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        /* Check detection of argument 2   */
        TRSYLV_CALL(2, "N", "X", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        /* Check detection of argument 3   */
        sign = 0;   TRSYLV_CALL(3, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        sign = 2;   TRSYLV_CALL(3, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        sign = -2;  TRSYLV_CALL(3, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        sign = 1;
        /* Check detection of argument 4  */
        M = -1; TRSYLV_CALL(4, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0;
        /* Check detection of argument 5  */
        N = -1; TRSYLV_CALL(5, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 0;
        /* Check detection of argument 7  */
        LDA = 0; TRSYLV_CALL(7, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 2; LDA = 1; TRSYLV_CALL(7, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0;
        /* Check detection of argument 9  */
        LDB = 0; TRSYLV_CALL(9, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 2; LDB = 1; TRSYLV_CALL(9, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 0;
        /* Check detection of argument 11  */
        LDX = 0; TRSYLV_CALL(11, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 2; LDA = 2;  LDX = 1; TRSYLV_CALL(11, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0; LDA = 1;  LDX = 1;

        if ( trsylv_functions[i].chk_l3 ) {
            mepack_double_trsylv_blocksize_mb_set(1);
            TRSYLV_CALL(30, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_trsylv_blocksize_mb_set(2);
            mepack_double_trsylv_blocksize_nb_set(1);
            TRSYLV_CALL(31, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_trsylv_blocksize_nb_set(2);
        }
        if ( trsylv_functions[i].chk_l3_2s ) {
            mepack_double_trsylv_blocksize_mb_set(4);
            mepack_double_trsylv_blocksize_mb_set(4);
            mepack_double_trsylv_blocksize_2stage_set(2);
            TRSYLV_CALL(33, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_trsylv_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

