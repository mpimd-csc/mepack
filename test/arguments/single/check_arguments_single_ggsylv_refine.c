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

typedef void (*tgsylv_function) (
            const char * TRANSA, const char * TRANSB, const char *GUESS, float *SGN, Int * M , Int* N,  float * A, Int * LDA,
            float *B, Int *LDB, float *C, Int * LDC, float * D, Int * LDD, float *X, Int *LDX, float *Y, Int *LDY,
            float *AS, Int * LDAS, float * BS, Int *LDBS, float * CS, Int * LDCS, float * DS, Int * LDDS, float * Q, Int *LDQ, float * Z, Int * LDZ,
            float *U, Int *LDU, float * V, Int * LDV,  Int * MAXIT, float *TAU, float *CONVLOG, float * WORK, Int * LWORK, Int *INFO,
            FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);


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

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,A34,A35,A36,A37,A38,A39,A40)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
    tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,A34,A35,A36,A37,A38,A39,A40, 1, 1, 1); \
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
    float WORK = 0 , sign = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, LDX=0, INFO=0;
    Int LDY = 0;
    Int LDC = 0, LDD = 0;
    Int LDQA=0, LDZA =0 , LDQB =0 , LDZB = 0 ;
    Int LDAS=0, LDBS=0, LDCS = 0, LDDS= 0;
    Int LDWORK = 1024*1024;
    float TAU =0;
    Int MAXIT = 0;
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
        { &(FC_GLOBAL_(sla_ggsylv_refine,SLA_GGSYLV_REFINE)), "SLA_GGSYLV_REFINE", 1, 0,0 },
        { NULL, NULL, 0, 0, 0}
    };


    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( tgsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", tgsylv_functions[i].name);
        tgsylv = tgsylv_functions[i].fn;

        /* Argument 1   */
        TGSYLV_CALL(1, "X", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 2   */
        TGSYLV_CALL(2, "N", "X", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 3   */
        TGSYLV_CALL(3, "N", "N", "G",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 4   */
        sign =-2;
        TGSYLV_CALL(4, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 5   */
        sign = 1;
        M = -1;
        TGSYLV_CALL(5, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 6  */
        sign = 1;
        M = 1;
        N = -1;
        TGSYLV_CALL(6, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 8   */
        sign = 1;
        N = 1;
        LDA = 0;
        TGSYLV_CALL(8, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 10   */
        sign = 1;
        N = 1;
        LDA = 1;
        LDB = 0;
        TGSYLV_CALL(10, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 12   */
        sign = 1;
        N = 1;
        LDC = 0;
        LDB = 1;
        TGSYLV_CALL(12, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 14   */
        sign = 1;
        N = 1;
        LDD = 0;
        LDC = 1;
        TGSYLV_CALL(14, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 16   */
        sign = 1;
        N = 1;
        LDX = 0;
        LDD = 1;
        TGSYLV_CALL(16, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 18   */
        sign = 1;
        N = 1;
        LDY = 0;
        LDX = 1;
        TGSYLV_CALL(18, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);



        /* Argument 20   */
        sign = 1;
        N = 1;
        LDAS = 0;
        LDY = 1;
        TGSYLV_CALL(20, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 22  */
        sign = 1;
        N = 1;
        LDBS = 0;
        LDAS = 1;
        TGSYLV_CALL(22, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 24   */
        sign = 1;
        N = 1;
        LDCS = 0;
        LDBS = 1;
        TGSYLV_CALL(24, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 26  */
        sign = 1;
        N = 1;
        LDDS = 0;
        LDCS = 1;
        TGSYLV_CALL(26, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 28   */
        sign = 1;
        N = 1;
        LDQA = 0;
        LDDS = 1;
        TGSYLV_CALL(28, "N", "N",  "N", &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 30   */
        sign = 1;
        N = 1;
        LDZA = 0;
        LDQA = 1;
        TGSYLV_CALL(30, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 32   */
        sign = 1;
        N = 1;
        LDQB = 0;
        LDZA = 1;
        TGSYLV_CALL(32, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 34   */
        sign = 1;
        N = 1;
        LDZB = 0;
        LDQB = 1;
        TGSYLV_CALL(34, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 35   */
        sign = 1;
        N = 1;
        LDZB = 1;
        MAXIT = 1;
        TGSYLV_CALL(35, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

       /* Argument 36   */
        sign = 1;
        N = 1;
        LDZB = 1;
        MAXIT = 40;
        TAU = -1.0;
        TGSYLV_CALL(36, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

       /* Argument 38   */
        sign = 1;
        N = 1;
        LDZB = 1;
        MAXIT = 41;
        TAU = 1.0;
        LDWORK = 0;
        TGSYLV_CALL(38, "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        printf("\n");
        LDWORK = 1024*1024;
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

