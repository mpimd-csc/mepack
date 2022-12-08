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

typedef void (*trlyap_function) (
            const char * TRANSA, const char *GUESS, Int * M , double * A, Int * LDA,
            double *X, Int *LDX, double *Y, Int *LDY,
            double *AS, Int * LDAS, double * Q, Int *LDQ,
            Int * MAXIT, double *TAU, double *CONVLOG, double * WORK, Int * LWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


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

#define TRLYAP_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15, A16, A17, A18, A19)  do { cname = trlyap_functions[i].name; global_state =0 ;  \
    ign_wrng = trlyap_functions[i].ign_wrng_name, count_all++;  \
    trlyap(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,1,1); \
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
    Int M = 0, LDA=0,INFO=0;
    Int LDQA=0;
    Int LDAS=0;
    Int LDX= 0;
    Int LDY=0;
    Int LDWORK = 1024*1024;
    double TAU =0;
    Int MAXIT = 0;
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
        { &(FC_GLOBAL_(dla_gelyap_refine,DLA_GELYAP_REFINE)), "DLA_GELYAP_REFINE", 1, 0,0 },
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

        /* Argument 1 - TRANS  */
        TRLYAP_CALL(1, "X", "N", &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        /* Argument 2 - Guess  */
        TRLYAP_CALL(2, "N", "X", &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 3 - M  */
        M = -1;
        TRLYAP_CALL(3, "N", "N", &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 5 - LDA  */
        M = 1;
        LDA = 0;
        TRLYAP_CALL(5, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);


        /* Argument 7 - LDX   */
        M = 1;
        LDA = 1;
        LDX = 0;
        TRLYAP_CALL(7, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);
        /* Argument 9 - LDY   */
        M = 1;
        LDA = 1;
        LDX = 1;
        LDY = 0;
        TRLYAP_CALL(9, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 11 - LDAS  */
        M = 1;
        LDX = 1 ;
        LDY = 1 ;
        LDAS = 0;
        TRLYAP_CALL(11, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 13 - LDQ  */
        M = 1;
        LDAS = 1 ;
        LDQA = 0 ;
        TRLYAP_CALL(13, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);


        /* Argument 14 - MAXIT  */
        M = 1;
        LDQA = 1 ;
        MAXIT = 1;

        TRLYAP_CALL(14, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 15 - TAU   */
        M = 1;
        MAXIT = 2 ;
        TAU = -1;
        TRLYAP_CALL(15, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        /* Argument 17  - LDWORK  */
        M = 1;
        TAU = 1 ;
        LDWORK = 0;
        TRLYAP_CALL(17, "N","N",  &M, NULL, &LDA, NULL, &LDX, NULL, &LDY,
                NULL, &LDAS,
                NULL, &LDQA, &MAXIT, &TAU, NULL, &WORK, &LDWORK, &INFO);

        LDWORK = 1000;




        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

