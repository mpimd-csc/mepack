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
        const char * TRANSA, const char *GUESS, Int * M , float * A, Int * LDA,
        float *B, Int *LDB, float *X, Int *LDX, float *Y, Int *LDY,
        float *AS, Int * LDAS, float * BS, Int *LDBS, float * Q, Int *LDQ, float * Z, Int * LDZ,
        Int * MAXIT, float *TAU, float *CONVLOG,  float * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


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

#define TGSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25)  do { cname = tgsylv_functions[i].name; global_state =0 ;  \
    ign_wrng = tgsylv_functions[i].ign_wrng_name, count_all++;  \
    tgsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24,A25,1,1); \
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
    float WORK = 0;
    Int M = 0, LDA=0, LDB=0, INFO=0;
    Int LDQA=0, LDZA =0 ;
    Int LDAS=0, LDBS=0;
    Int LDX= 0;
    Int LDY= 0;
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
        { &(FC_GLOBAL_(sla_gglyap_refine,SLA_GGLYAP_REFINE)), "SLA_GGLYAP_REFINE", 1, 0,0 },
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

        /* Argument 1   */
        TGSYLV_CALL(1, "X", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 2   */
        TGSYLV_CALL(2, "N", "G", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);


        /* Argument 3   */
        M = -1;
        TGSYLV_CALL(3, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 5   */
        M = 1;
        LDA = 0;
        TGSYLV_CALL(5, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 7   */
        M = 1;
        LDA = 1 ;
        LDB = 0;
        TGSYLV_CALL(7, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 9   */
        M = 1;
        LDB = 1 ;
        LDX = 0;
        TGSYLV_CALL(9, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 11  */
        M = 1;
        LDX = 1 ;
        LDY = 0;
        TGSYLV_CALL(11, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);



        /* Argument 13   */
        M = 1;
        LDY = 1 ;
        LDAS = 0;
        TGSYLV_CALL(13, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 15  */
        M = 1;
        LDAS = 1 ;
        LDBS = 0;
        TGSYLV_CALL(15, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 17   */
        M = 1;
        LDBS = 1 ;
        LDQA = 0;
        TGSYLV_CALL(17, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 19   */
        M = 1;
        LDQA = 1 ;
        LDZA = 0;
        TGSYLV_CALL(19, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 20   */
        M = 1;
        LDQA = 1 ;
        LDZA = 1;
        MAXIT = 1;
        TGSYLV_CALL(20, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 21   */
        M = 1;
        MAXIT = 2 ;
        TAU = -1;
        TGSYLV_CALL(21, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 23   */
        M = 1;
        TAU = 1 ;
        LDWORK = 0;
        TGSYLV_CALL(23, "N", "N", &M, NULL, &LDA, NULL, &LDB, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS, NULL, &LDBS,
                NULL, &LDQA, NULL, &LDZA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        LDWORK = 1000*1024;




        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

