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
        double *sgn1, double *sgn2,  Int *M, Int *N, double * A, Int *LDA,double * B, Int * LDB,
        double *C, Int *ldc,double* D, Int *LDd,
        double *QA, Int *LDQA, double *ZA, Int *LDZA, double *QB, Int *LDQB, double *ZB, Int *LDZB,
        double *E, Int *LDE, double *F, Int *LDF, double * SCALE, double *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2, FORTRAN_LEN_T l3, FORTRAN_LEN_T l4);


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

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
    tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,1,1,1,1); \
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
    double WORK = 0 , sign1 = 0, sign2=0, SCALE = 0;
    Int M = 0, N =0 , LDA=0, LDB=0, LDE=0, LDF=0,  INFO=0;
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
        int dual;
    };

    struct function_name tgsylv_functions[] = {
        { &(FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)), "DLA_GGCSYLV", 1, MEPACK_FRONTEND_SOLVER_LEVEL3,0, 0 },
        { &(FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)), "DLA_GGCSYLV", 1, MEPACK_FRONTEND_SOLVER_LEVEL2,0,0 },
        { &(FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)), "DLA_GGCSYLV", 1, MEPACK_FRONTEND_SOLVER_DAG,0,0 },
        { &(FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)), "DLA_GGCSYLV", 1, MEPACK_FRONTEND_SOLVER_2STAGE,0,0 },
        { &(FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)), "DLA_GGCSYLV", 1, MEPACK_FRONTEND_SOLVER_RECURSIVE,0,0 },
        { &(FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)), "DLA_GGCSYLV_DUAL", 1, MEPACK_FRONTEND_SOLVER_LEVEL3,0, 1 },
        { &(FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)), "DLA_GGCSYLV_DUAL", 1, MEPACK_FRONTEND_SOLVER_LEVEL2,0, 1 },
        { &(FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)), "DLA_GGCSYLV_DUAL", 1, MEPACK_FRONTEND_SOLVER_DAG,0,1 },
        { &(FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)), "DLA_GGCSYLV_DUAL", 1, MEPACK_FRONTEND_SOLVER_2STAGE,0,1 },
        { &(FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)), "DLA_GGCSYLV_DUAL", 1, MEPACK_FRONTEND_SOLVER_RECURSIVE,0,1 },
        { NULL, NULL, 0, 0, 0, 0 },
    };


    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( tgsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", tgsylv_functions[i].name);
        tgsylv = tgsylv_functions[i].fn;
        if ( tgsylv_functions[i].dual ) {
            mepack_tgcsylv_dual_frontend_solver_set(tgsylv_functions[i].isolver);
        } else {
            mepack_tgcsylv_frontend_solver_set(tgsylv_functions[i].isolver);
        }
        /* Argument 1   */
        TGSYLV_CALL(1, "X", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        /* Argument 2   */
        TGSYLV_CALL(2, "N", "X", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        /* Argument 3   */
        TGSYLV_CALL(3, "N", "N", "X", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        /* Argument 4   */
        TGSYLV_CALL(4, "N", "N", "N", "X",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        /* Argument 5  */
        sign1 = 0;
        TGSYLV_CALL(5, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        sign1 = -2;
        TGSYLV_CALL(5, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
         sign1 = 2;
        TGSYLV_CALL(5, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        sign1 = 1;

        /* Argument 6  */
        sign2 = 0;
        TGSYLV_CALL(6, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        sign2 = -2;
        TGSYLV_CALL(6, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
         sign2 = 2;
        TGSYLV_CALL(6, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        sign2 = 1;

        /* Argument 7   */
        M = -1;
        TGSYLV_CALL(7, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        M = 0;
        /* Argument 8    */
        N = -1;
        TGSYLV_CALL(8, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        N = 0;

        /* Argument 10  */
        LDA = 0;
        TGSYLV_CALL(10, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDA = 1; M = 2;
        TGSYLV_CALL(10, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDA = 5;

        /* Argument 12 */
        LDB = 0;
        TGSYLV_CALL(12, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDB = 1; N = 2;
        TGSYLV_CALL(12, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDB = 5;

        /* Argument 14  */
        LDC = 0;
        TGSYLV_CALL(14, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDC = 1; M = 2;
        TGSYLV_CALL(14, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDC = 5;

        /* Argument 16 */
        LDD = 0;
        TGSYLV_CALL(16, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDD = 1; N = 2;
        TGSYLV_CALL(16, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDD = 5;

        /* Argument 18  */
        LDQA = 0;
        TGSYLV_CALL(18, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDQA = 1; M = 2;
        TGSYLV_CALL(18, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDQA = 5;

        /* Argument 20  */
        LDZA = 0;
        TGSYLV_CALL(20, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDZA = 1; M = 2;
        TGSYLV_CALL(20, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDZA = 5;

        /* Argument 22 */
        LDQB = 0;
        TGSYLV_CALL(22, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDQB = 1; N = 2;
        TGSYLV_CALL(22, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDQB = 5;

        /* Argument 24 */
        LDZB = 0;
        TGSYLV_CALL(24, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDZB = 1; N = 2;
        TGSYLV_CALL(24, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDZB = 5;
        /* Argument 26 */
        LDE = 0;
        TGSYLV_CALL(26, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDE = 1; M = 2;
        TGSYLV_CALL(26, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDE = 5;
        /* Argument 28 */
        LDF = 0;
        TGSYLV_CALL(28, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);

        LDF = 1; M = 2;
        TGSYLV_CALL(28, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDF = 5;


        /* Argument 31 */
        LDWORK = 0;
        TGSYLV_CALL(31, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
        LDWORK = 1024*1024;

        if ( tgsylv_functions[i].chk_l3 ) {
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_mb_set(1); else mepack_double_tgcsylv_blocksize_mb_set(1);

            TGSYLV_CALL(40, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                    NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_mb_set(2); else mepack_double_tgcsylv_blocksize_mb_set(2);

            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_nb_set(1); else mepack_double_tgcsylv_blocksize_nb_set(1);
            TGSYLV_CALL(41, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                    NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_nb_set(2); else mepack_double_tgcsylv_blocksize_nb_set(2);

        }
        if ( tgsylv_functions[i].isolver == MEPACK_FRONTEND_SOLVER_2STAGE ) {
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_nb_set(4); else mepack_double_tgcsylv_blocksize_nb_set(4);
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_mb_set(4); else mepack_double_tgcsylv_blocksize_mb_set(4);
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_2stage_set(2); else mepack_double_tgcsylv_blocksize_2stage_set(2);
            TGSYLV_CALL(43, "N", "N", "N", "N",  &sign1, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,
                    NULL, &LDQA, NULL, &LDZA, NULL, &LDQB, NULL, &LDZB, NULL, &LDE, NULL, &LDF, &SCALE, &WORK, &LDWORK, &INFO);
            if ( tgsylv_functions[i].dual ) mepack_double_tgcsylv_dual_blocksize_2stage_set(512); else mepack_double_tgcsylv_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

