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
            const char * TRANSA, const char * TRANSB, const char *GUESS, double *SGN1, double *SGN2, Int * M , Int* N,  double * A, Int * LDA,
            double *B, Int *LDB, double *C, Int * LDC, double * D, Int * LDD, double *R, Int *LDR, double *L, Int *LDL, double *E, Int *LDE, double *F, Int *LDF,
            double *AS, Int * LDAS, double * BS, Int *LDBS, double * CS, Int * LDCS, double * DS, Int * LDDS, double * Q, Int *LDQ, double * Z, Int * LDZ,
            double *U, Int *LDU, double * V, Int * LDV,  Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO,
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

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,A34,A35,A36,A37,A38,A39,A40,A41,A42,A43,A44,A45)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
    tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,A34,A35,A36,A37,A38,A39,A40,A41,A42,A43,A44,A45,1,1,1); \
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
    double WORK = 0 , sign1 = 0, sign2 = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, INFO=0;
    Int LDC = 0, LDD = 0;
    Int LDQA=0, LDZA =0 , LDQB =0 , LDZB = 0 ;
    Int LDAS=0, LDBS=0, LDCS = 0, LDDS= 0;
    Int LDE=0, LDF = 0;
    Int LDR=0, LDL = 0;
    Int LDWORK = 1024*1024;
    double TAU =0;
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
        int dual;
    };

    struct function_name tgsylv_functions[] = {
        { &(FC_GLOBAL_(dla_ggcsylv_refine,DLA_GGCSYLV_REFINE)), "DLA_GGCSYLV_REFINE", 1, 0,0,0 },
        { &(FC_GLOBAL_(dla_ggcsylv_dual_refine,DLA_GGCSYLV_DUAL_REFINE)), "DLA_GGCSYLV_DUAL_REFINE", 1, 0,0,1 },
        { NULL, NULL, 0,0, 0, 0}
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
        TGSYLV_CALL(1, "X", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 2   */
        TGSYLV_CALL(2, "N", "X", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 3 "GUESS"*/
        TGSYLV_CALL(3, "N", "N", "G",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);


        /* Argument 4   */
        sign1 = -2;
        TGSYLV_CALL(4, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 5   */
        sign1 = -1;
        sign2 = -3;
        TGSYLV_CALL(5, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);


        /* Argument 6   */
        sign2 = 1;
        M = -1;
        TGSYLV_CALL(6, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 7  */
        M = 1;
        N = -1;
        TGSYLV_CALL(7, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 9   */
        N = 1;
        LDA = 0;
        TGSYLV_CALL(9, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 11   */
        N = 1;
        LDA = 1;
        LDB = 0;
        TGSYLV_CALL(11, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 13   */
        N = 1;
        LDC = 0;
        LDB = 1;
        TGSYLV_CALL(13, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 15   */
        N = 1;
        LDD = 0;
        LDC = 1;
        TGSYLV_CALL(15, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 17   */
        N = 1;
        LDR = 0;
        LDD = 1;
        TGSYLV_CALL(17, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 19   */
        N = 1;
        LDL = 0;
        LDR = 1;
        TGSYLV_CALL(19, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);



        /* Argument 21   */
        N = 1;
        LDE = 0;
        LDL = 1;
        TGSYLV_CALL(21, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 23   */
        N = 1;
        LDE = 1;
        LDF = 0;
        TGSYLV_CALL(23, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);



        /* Argument 25   */
        N = 1;
        LDE = 1;
        LDF = 1;
        LDAS = 0;
        TGSYLV_CALL(25, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 27  */
        N = 1;
        LDBS = 0;
        LDAS = 1;
        TGSYLV_CALL(27, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 29   */
        N = 1;
        LDCS = 0;
        LDBS = 1;
        TGSYLV_CALL(29, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 31  */
        N = 1;
        LDDS = 0;
        LDCS = 1;
        TGSYLV_CALL(31, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 33   */
        N = 1;
        LDQA = 0;
        LDDS = 1;
        TGSYLV_CALL(33, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 35   */
        N = 1;
        LDZA = 0;
        LDQA = 1;
        TGSYLV_CALL(35, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 37   */
        N = 1;
        LDQB = 0;
        LDZA = 1;
        TGSYLV_CALL(37, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 39   */
        N = 1;
        LDZB = 0;
        LDQB = 1;
        TGSYLV_CALL(39, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 40   */
        N = 1;
        LDZB = 1;
        MAXIT = 1;
        TGSYLV_CALL(40, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

       /* Argument 41  */
        N = 1;
        LDZB = 1;
        MAXIT = 2;
        TAU = -1.0;
        TGSYLV_CALL(41, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

       /* Argument 43   */
        N = 1;
        LDZB = 1;
        MAXIT = 2;
        TAU = 1.0;
        LDWORK = 0;
        TGSYLV_CALL(43, "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, NULL, &LDR, NULL, &LDL, NULL, &LDE , NULL, &LDF,
                NULL, &LDAS, NULL, &LDBS, NULL, &LDCS, NULL, &LDDS,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

