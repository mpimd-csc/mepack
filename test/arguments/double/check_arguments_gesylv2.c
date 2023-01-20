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
typedef void (*trsylv_function) (const char *FACTA, const char* FACTB, const char *TRANSA, const char*TRANSB, double *sgn,
        Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB,
        double *QA, Int *LDQA,double * QB, Int * LDQB,
        double *X, Int *LDX, double * SCALE,
        double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);


Int global_state = 0;
char *cname = NULL;
int ign_wrng;

void err_handler(const char * str, Int Info)
{
    if ( strcmp(str,cname) == 0 || ign_wrng )
        global_state = Info;
    else {
        printf("-- Wrong Name in error handler ( %s) .\n", str);
        global_state = -Info;
    }
}

#define TRSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21)  do { cname = trsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = trsylv_functions[i].ign_wrng_name, count_all++;  \
    trsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21, 1,1,1,1); \
    if ( global_state == D ) {\
        ok &= 1; count_ok ++;  \
        printf("-- ARGUMENT %d OK.\n",D); \
    } else { \
        printf("-- ARGUMENT %d FAILED. INFO = %d\n",D,(int) -global_state); \
        ok &=0; count_failed++;   \
    } \
} while(0);

int main(int argc, char **argv)
{
    double WORK = 0 , sign = 0, SCALE = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, LDX=0, INFO=0;
    Int LDQA = 0, LDQB =1, LDWORK = 1;
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
        int isolver;
        int chk_isolve;
        int ign_wrng_name;
    };

    struct function_name trsylv_functions[] = {
        { &(FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)), "DLA_GESYLV2", 1, MEPACK_FRONTEND_SOLVER_LEVEL3, 1, 0 },
        { &(FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)), "DLA_GESYLV2", 1, MEPACK_FRONTEND_SOLVER_LEVEL2, 1, 0 },
        { &(FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)), "DLA_GESYLV2", 1, MEPACK_FRONTEND_SOLVER_DAG, 1, 0 },
        { &(FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)), "DLA_GESYLV2", 1, MEPACK_FRONTEND_SOLVER_2STAGE, 1, 0 },
        { &(FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)), "DLA_GESYLV2", 1, MEPACK_FRONTEND_SOLVER_RECURSIVE, 1, 0 },
        { NULL, NULL, 0,0, 0, 0}
    };


    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( trsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", trsylv_functions[i].name);
        mepack_trsylv2_frontend_solver_set(trsylv_functions[i].isolver);
        trsylv = trsylv_functions[i].fn;

        /* Check detection of argument 1   */
        TRSYLV_CALL(1, "X", "N", "X", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        /* Check detection of argument 2 */
        TRSYLV_CALL(2, "N", "X", "X", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        /* Check detection of argument 3 */
        TRSYLV_CALL(3, "N", "N", "X", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        /* Check detection of argument 4 */
        TRSYLV_CALL(4, "N", "N", "N", "X", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        /* Check detection of argument 5 */
        sign = 0; TRSYLV_CALL(5, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        sign = 2; TRSYLV_CALL(5, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        sign = -2; TRSYLV_CALL(5, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        sign = 1;
        /* Check detection of argument 6 */
        M=-1; TRSYLV_CALL(6, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        M = 0;
        /* Check detection of argument 7 */
        N=-1; TRSYLV_CALL(7, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        N = 0;

        /* Check detection of argument 9 */
        LDA=0; TRSYLV_CALL(9, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        M=2; LDA=1; TRSYLV_CALL(9, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        LDA = 5;

        /* Check detection of argument 11 */
        LDB=0; TRSYLV_CALL(11, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        N=2; LDB=1; TRSYLV_CALL(11, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        LDB = 5;

        /* Check detection of argument 13 */
        LDQA=0; TRSYLV_CALL(13, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        M=2; LDQA=1; TRSYLV_CALL(13, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        LDQA = 5;

        /* Check detection of argument 15 */
        LDQB=0; TRSYLV_CALL(15, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        N=2; LDQB=1; TRSYLV_CALL(15, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        LDQB = 5;

        /* Check detection of argument 17 */
        LDX=0; TRSYLV_CALL(17, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        M=2; LDX=1; TRSYLV_CALL(17, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        LDX = 5;
        /* Check detection of argument 20   */
        M=N=2; LDWORK=0; TRSYLV_CALL(20, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
        LDWORK = 2000*2000;


        if ( trsylv_functions[i].chk_l3 ) {
            mepack_double_trsylv2_blocksize_mb_set(1);
            TRSYLV_CALL(30, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
            mepack_double_trsylv2_blocksize_mb_set(2);

            mepack_double_trsylv2_blocksize_nb_set(1);
            TRSYLV_CALL(31, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
            mepack_double_trsylv2_blocksize_nb_set(2);
        }
        if ( trsylv_functions[i].isolver == MEPACK_FRONTEND_SOLVER_2STAGE ) {
            mepack_double_trsylv2_blocksize_mb_set(4);
            mepack_double_trsylv2_blocksize_nb_set(4);
            mepack_double_trsylv2_blocksize_2stage_set(2);
            TRSYLV_CALL(33, "N", "N", "N", "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDQA, NULL, &LDQB, NULL, &LDX, &SCALE, &WORK, &LDWORK,  &INFO);
            mepack_double_trsylv2_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

