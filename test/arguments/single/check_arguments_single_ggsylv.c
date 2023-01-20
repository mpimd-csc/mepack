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

typedef void (*tgsylv_function) (const char * FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB,
        float *sgn,  Int *M, Int *N, float * A, Int *LDA,float * B, Int * LDB,
        float *C, Int *ldc,float* D, Int *LDd,
        float *QA, Int *LDQA, float *ZA, Int *LDZA, float *QB, Int *LDQB, float *ZB, Int *LDZB,
        float *X, Int *LDX, float * SCALE, float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);


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

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
    tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29, 1, 1, 1, 1); \
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
    float WORK = 0 , sign = 0, SCALE = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, LDX=0, INFO=0;
    Int LDC = 0, LDD = 0;
    Int LDQA=0, LDZA =0 , LDQB =0 , LDZB = 0 ;
    Int LDWORK = 1024*1024;
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
        int isolver;
        int ign_wrng_name;
    };

    struct function_name tgsylv_functions[] = {
        { &(FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)), "SLA_GGSYLV", 1, MEPACK_FRONTEND_SOLVER_LEVEL3,0 },
        { &(FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)), "SLA_GGSYLV", 1, MEPACK_FRONTEND_SOLVER_LEVEL2,0 },
        { &(FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)), "SLA_GGSYLV", 1, MEPACK_FRONTEND_SOLVER_DAG,0 },
        { &(FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)), "SLA_GGSYLV", 1, MEPACK_FRONTEND_SOLVER_2STAGE,0 },
        { &(FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)), "SLA_GGSYLV", 1, MEPACK_FRONTEND_SOLVER_RECURSIVE,0 },
        { &(FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)), "SLA_GGSYLV", 1, MEPACK_FRONTEND_SOLVER_GARDINER_LAUB,0 },
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

        mepack_tgsylv_frontend_solver_set(tgsylv_functions[i].isolver);

        /* Argument 1   */
        TGSYLV_CALL(1, "X", "N", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        /* Argument 2   */
        TGSYLV_CALL(2, "F", "X", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        /* Argument 3   */
        TGSYLV_CALL(3, "F", "F", "X", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        /* Argument 4   */
        TGSYLV_CALL(4, "F", "F", "N", "X",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        /* Argument 5  */
        sign = 0;
        TGSYLV_CALL(5, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        sign = -2;
        TGSYLV_CALL(5, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        sign = 2;
        TGSYLV_CALL(5, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        sign = 1;

        /* Argument 6   */
        M = -1;
        TGSYLV_CALL(6, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M = 0;
        /* Argument 7    */
        N = -1;
        TGSYLV_CALL(7, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        N = 0;

        /* Argument 9  */
        LDA = 0;
        TGSYLV_CALL(9, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDA = 1; M = 2;
        TGSYLV_CALL(9, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDA = 5;

        /* Argument 11 */
        LDB = 0;
        TGSYLV_CALL(11, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDB = 1; N = 2;
        TGSYLV_CALL(11, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDB = 5;

        /* Argument 13  */
        LDC = 0;
        TGSYLV_CALL(13, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDC = 1; M = 2;
        TGSYLV_CALL(13, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDC = 5;

        /* Argument 15 */
        LDD = 0;
        TGSYLV_CALL(15, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDD = 1; N = 2;
        TGSYLV_CALL(15, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDD = 5;

        /* Argument 17  */
        LDQA = 0;
        TGSYLV_CALL(17, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDQA = 1; M = 2;
        TGSYLV_CALL(17, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDQA = 5;

        /* Argument 19  */
        LDZA = 0;
        TGSYLV_CALL(19, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDZA = 1; M = 2;
        TGSYLV_CALL(19, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDZA = 5;

        /* Argument 21 */
        LDQB = 0;
        TGSYLV_CALL(21, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDQB = 1; N = 2;
        TGSYLV_CALL(21, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDQB = 5;

        /* Argument 23 */
        LDZB = 0;
        TGSYLV_CALL(23, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDZB = 1; N = 2;
        TGSYLV_CALL(23, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDZB = 5;

        /* Argument 25 */
        LDX = 0;
        TGSYLV_CALL(25, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        LDX = 1; M = 2;
        TGSYLV_CALL(25, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDX = 5;

        /* Argument 28 */
        LDWORK = 0;
        TGSYLV_CALL(28, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDWORK = 1024*1024;

        if ( tgsylv_functions[i].chk_l3 ) {
            mepack_single_tgsylv_blocksize_mb_set(1);
            TGSYLV_CALL(30, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                    NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_single_tgsylv_blocksize_mb_set(2);
            mepack_single_tgsylv_blocksize_nb_set(1);
            TGSYLV_CALL(31, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                    NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_single_tgsylv_blocksize_nb_set(2);
            /* mepack_tgsylv_isolver_set(0); */
            /* TGSYLV_CALL(32, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, */
            /*         NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO); */
            /*  */
            /* mepack_tgsylv_isolver_set(7); */
            /* TGSYLV_CALL(32, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD, */
            /*         NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO); */
            /*  */
            /* mepack_tgsylv_isolver_set(1); */
        }
        if ( tgsylv_functions[i].isolver == MEPACK_FRONTEND_SOLVER_2STAGE ) {
            mepack_single_tgsylv_blocksize_mb_set(4);
            mepack_single_tgsylv_blocksize_nb_set(4);
            mepack_single_tgsylv_blocksize_2stage_set(2);
            TGSYLV_CALL(33, "F", "F", "N", "N",  &sign, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                    NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_single_tgsylv_blocksize_2stage_set(512);
        }


        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

