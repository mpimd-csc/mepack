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

typedef void (*tgsylv_function) (const char *TRANSA, const char*TRANSB, double *sgn,  Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB, double *C, Int *ldc,double* D, Int *LDC, double *X, Int *LDX, double * SCALE, double *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


Int global_state = 0;
char *cname = NULL;
int ign_wrng;

void err_handler(const char * str, Int Info)
{
    if ( strcmp(str,cname) == 0 || ign_wrng )
        global_state = Info;
    else {
        printf("-- Wrong Name in error handler. str = %s  cname = %s\n", str, cname );
        global_state = -Info;
    }
}

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
            ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
            tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,1, 1); \
            *(A18) = 0; \
            if ( global_state == D ) {\
                ok &= 1; count_ok ++;  \
                printf("-- ARGUMENT %d OK.\n",D); \
            } else { \
                printf("-- ARGUMENT %d FAILED (%d).\n",D, (int) global_state); \
                ok &=0; count_failed++;   \
            } \
        } while(0);

int main(int argc, char **argv)
{
    double WORK = 0 , sign = 0, SCALE = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, LDX=0, INFO=0;
    Int LDC = 0, LDD = 0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    tgsylv_function tgsylv;
    struct function_name {
        tgsylv_function fn;
        char * name;
        int chk_l3;
        int chk_l3_2s;
        int ign_wrng_name;
    };

    struct function_name tgsylv_functions[] = {
        { &(FC_GLOBAL_(dla_tgsylv_l3,DLA_TGSYLV_L3)), "DLA_TGSYLV_L3", 1, 0,0 },
        { &(FC_GLOBAL_(dla_tgsylv_l3_rowwise,DLA_TGSYLV_L3_ROWWISE)), "DLA_TGSYLV_L3_ROWWISE", 1, 0,0  },
        { &(FC_GLOBAL_(dla_tgsylv_l3_colwise,DLA_TGSYLV_L3_COLWISE)), "DLA_TGSYLV_L3_COLWISE", 1, 0,0  },
        { &(FC_GLOBAL_(dla_tgsylv_l2,DLA_TGSYLV_L2)), "DLA_TGSYLV_L2", 0, 0 ,0 },
        { &(FC_GLOBAL_(dla_tgsylv_l2_rowwise,DLA_TGSYLV_L2_ROWWISE)), "DLA_TGSYLV_L2_ROWWISE", 0, 0 ,0 },
        { &(FC_GLOBAL_(dla_tgsylv_l2_colwise,DLA_TGSYLV_L2_COLWISE)), "DLA_TGSYLV_L2_COLWISE", 0, 0 ,0 },
        { &(FC_GLOBAL_(dla_tgsylv_l2_reorder,DLA_TGSYLV_L2_REORDER)), "DLA_TGSYLV_L2_REORDER", 0, 0,0  },
        { &(FC_GLOBAL_(dla_tgsylv_l2_local_copy,DLA_TGSYLV_L2_LOCAL_COPY)), "DLA_TGSYLV_L2_LOCAL_COPY" ,0 ,0,0  },
        { &(FC_GLOBAL_(dla_tgsylv_l2_local_copy_32,DLA_TGSYLV_L2_LOCAL_COPY_32)), "DLA_TGSYLV_L2_LOCAL_COPY_32", 0,0 , 1},
        { &(FC_GLOBAL_(dla_tgsylv_l2_local_copy_64,DLA_TGSYLV_L2_LOCAL_COPY_64)), "DLA_TGSYLV_L2_LOCAL_COPY_64", 0,0 , 1},
        { &(FC_GLOBAL_(dla_tgsylv_l2_local_copy_96,DLA_TGSYLV_L2_LOCAL_COPY_96)), "DLA_TGSYLV_L2_LOCAL_COPY_96", 0,0 , 1},
        { &(FC_GLOBAL_(dla_tgsylv_l2_local_copy_128,DLA_TGSYLV_L2_LOCAL_COPY_128)), "DLA_TGSYLV_L2_LOCAL_COPY_128", 0,0,1 },
        { &(FC_GLOBAL_(dla_tgsylv_l2_unopt,DLA_TGSYLV_L2_UNOPT)), "DLA_TGSYLV_L2_UNOPT", 0,0, 0 },
        { &(FC_GLOBAL_(dla_tgsylv_dag,DLA_TGSYLV_DAG)), "DLA_TGSYLV_DAG", 1, 0, 0 },
        { &(FC_GLOBAL_(dla_tgsylv_l3_2s,DLA_TGSYLV_L3_2S)), "DLA_TGSYLV_L3_2S",1 , 1 ,0},
        { &(FC_GLOBAL_(dla_tgsylv_recursive,DLA_TGSYLV_RECURSIVE)), "DLA_TGSYLV_RECURSIVE", 0, 0 , 0},
        { NULL, NULL, 0,0, 0}
    };

    argc = argc;
    (void) argv;


    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( tgsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", tgsylv_functions[i].name);
        tgsylv = tgsylv_functions[i].fn;

        /* Check detection of argument 1   */
        TGSYLV_CALL(1, "X", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        /* Check detection of argument 2   */
        TGSYLV_CALL(2, "N", "X", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        /* Check detection of argument 3   */
        sign = 0;   TGSYLV_CALL(3, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        sign = 2;   TGSYLV_CALL(3, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        sign = -2;  TGSYLV_CALL(3, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        sign = 1;
        /* Check detection of argument 4  */
        M = -1; TGSYLV_CALL(4, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0;
        /* Check detection of argument 5  */
        N = -1; TGSYLV_CALL(5, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 0;
        /* Check detection of argument 7  */
        LDA = 0;        TGSYLV_CALL(7, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 2; LDA = 1; TGSYLV_CALL(7, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0; LDA = 5;
        /* Check detection of argument 9  */
        LDB = 0; TGSYLV_CALL(9, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 2; LDB = 1; TGSYLV_CALL(9, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 0; LDB = 5;

        /* Check detection of argument 11  */
        LDC = 0;        TGSYLV_CALL(11, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 2; LDC = 1; TGSYLV_CALL(11, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0; LDC = 5;

        /* Check detection of argument 13  */
        LDD = 0; TGSYLV_CALL(13, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 2; LDD = 1; TGSYLV_CALL(13, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        N = 0; LDD = 5;

        /* Check detection of argument 15  */
        LDX = 0; TGSYLV_CALL(15, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDC, NULL, &LDD, NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 2; LDA = 2;  LDX = 1; TGSYLV_CALL(15, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
        M = 0; LDA = 1;  LDX = 1;

        if ( tgsylv_functions[i].chk_l3 ) {
            mepack_double_tgsylv_blocksize_mb_set(1);
            TGSYLV_CALL(30, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_tgsylv_blocksize_mb_set(2);
            mepack_double_tgsylv_blocksize_nb_set(1);
            TGSYLV_CALL(31, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_tgsylv_blocksize_nb_set(2);
        }
        if ( tgsylv_functions[i].chk_l3_2s ) {
            mepack_double_tgsylv_blocksize_mb_set(4);
            mepack_double_tgsylv_blocksize_nb_set(4);
            mepack_double_tgsylv_blocksize_2stage_set(2);
            TGSYLV_CALL(33, "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, &SCALE, &WORK, &INFO);
            mepack_double_tgsylv_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

