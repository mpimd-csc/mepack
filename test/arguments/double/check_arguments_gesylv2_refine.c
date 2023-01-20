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

typedef void (*tgsylv_function) (
            const char * TRANSA, const char * TRANSB, const char *GUESS, double * SIGN, Int * M , Int *N, double * A, Int * LDA,
            double *B, Int *LDB, double *X, Int *LDX, double *Y, Int *LDY,
            double *AS, Int * LDAS, double * BS, Int *LDBS, double * Q, Int *LDQ, double * Z, Int * LDZ,
            Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3);


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

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22, A23, A24, A25, A26, A27, A28)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
    tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22, A23, A24, A25, A26, A27, A28, 1, 1, 1); \
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
    double WORK = 0;
    Int M = 0, N = 0, LDA=0, LDB=0, INFO=0;
    Int LDQA=0, LDZA =0 ;
    Int LDAS=0, LDBS=0;
    Int LDX= 0;
    Int LDY = 0;
    Int LDWORK = 1024*1024;
    double TAU =0, SGN = 1.0;
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
        { &(FC_GLOBAL_(dla_gesylv2_refine,DLA_GESYLV2_REFINE)), "DLA_GESYLV2_REFINE", 1, 0,0 },
        { NULL, NULL, 0,0,0}
    };

    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( tgsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", tgsylv_functions[i].name);
        tgsylv = tgsylv_functions[i].fn;

        /* Argument 1   TRANSA */
        TGSYLV_CALL(1, "X", "X", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        /* Argument 2   TRANSB */
        TGSYLV_CALL(2, "N", "X", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 3   GUESS */
        TGSYLV_CALL(3, "N", "N", "G", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);


        /* Argument 4   SIGN*/
        SGN = -10;
        TGSYLV_CALL(4, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        SGN = 1.0;


        /* Argument 5 - M   */
        M = -1;
        TGSYLV_CALL(5, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        M = 1;

        /* Argument 6  - N   */
        N = -1;
        TGSYLV_CALL(6, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        N = 0;

        /* Argument 8 - LDA   */
        LDA = 0;
        TGSYLV_CALL(8, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDA = 1;

        /* Argument 10 - LDB   */
        LDB = 0;
        TGSYLV_CALL(10, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDB = 1;

        /* Argument 12 - LDX   */
        LDX = 0;
        TGSYLV_CALL(12, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDX = 1;

        /* Argument 14 - LDY   */
        LDY = 0;
        TGSYLV_CALL(14, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDY = 1;



        /* Argument 16 - LDAS   */
        LDAS = 0;
        TGSYLV_CALL(16, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDAS = 1;

        /* Argument 18 - LDBS   */
        LDBS = 0;
        TGSYLV_CALL(18, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDBS = 1;

        /* Argument 20 - LDQA   */
        LDQA = 0;
        TGSYLV_CALL(20, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDQA = 1;

        /* Argument 22 - LDZA   */
        LDZA = 0;
        TGSYLV_CALL(22, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDZA = 1;

        /* Argument 23 - MAXIT   */
        MAXIT = 1;
        TGSYLV_CALL(23, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        MAXIT = 50;

        /* Argument 24 - TAU   */
        TAU = -1;
        TGSYLV_CALL(24, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        TAU = 1;

        /* Argument 26 - LDWORK   */
        LDWORK = 0;
        TGSYLV_CALL(26, "N", "N", "N", &SGN, &M, &N,  NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDWORK = 1024*1024;

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

