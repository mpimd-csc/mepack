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
typedef void (*trstein_function) (const char*FACT, const char *TRANSA, Int *M, float * A, Int *LDA, float *Q, Int *LDQ, float *X, Int *LDX, float * SCALE,
        float *WORK, Int *LDWORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


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

#define GESTEIN_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13)  do { cname = trstein_functions[i].name; global_state =0;  \
            ign_wrng = trstein_functions[i].ign_wrng_name, count_all++;  \
            trstein(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13, 1, 1); \
            if ( global_state == D ) {\
                ok &= 1; count_ok ++;  \
                printf("-- ARGUMENT %d OK.\n",D); \
            } else { \
                printf("-- ARGUMENT %d FAILED (info = %d).\n",D, (int) global_state); \
                ok &=0; count_failed++;   \
            } \
        } while(0);

int main(int argc, char **argv)
{
    float WORK = 0 , SCALE = 0;
    Int M = 0, LDA=0, LDX=0, INFO=0, LDQ = 0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    Int LDWORK;
    trstein_function trstein;
    struct function_name {
        trstein_function fn;
        char * name;
        int chk_l3;
        int isolver;
        int ign_wrng_name;
    };

    struct function_name trstein_functions[] = {
        { &(FC_GLOBAL_(sla_gestein,SLA_GESTEIN)), "SLA_GESTEIN", 1, MEPACK_FRONTEND_SOLVER_LEVEL3,0 },
        { &(FC_GLOBAL_(sla_gestein,SLA_GESTEIN)), "SLA_GESTEIN", 1, MEPACK_FRONTEND_SOLVER_LEVEL2,0 },
        { &(FC_GLOBAL_(sla_gestein,SLA_GESTEIN)), "SLA_GESTEIN", 1, MEPACK_FRONTEND_SOLVER_DAG,0 },
        { &(FC_GLOBAL_(sla_gestein,SLA_GESTEIN)), "SLA_GESTEIN", 1, MEPACK_FRONTEND_SOLVER_2STAGE,0 },
        { &(FC_GLOBAL_(sla_gestein,SLA_GESTEIN)), "SLA_GESTEIN", 1, MEPACK_FRONTEND_SOLVER_RECURSIVE,0 },
        { NULL, NULL, 0,0, 0}
    };


    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( trstein_functions[i].fn != NULL ) {
        printf("Check %s\n", trstein_functions[i].name);
        mepack_trstein_frontend_solver_set(trstein_functions[i].isolver);
        LDWORK = 4000;
        trstein = trstein_functions[i].fn;

        /* Check detection of argument 1   */
        GESTEIN_CALL(1, "B", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        /* Check detection of argument 2   */
        GESTEIN_CALL(2, "F", "X", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);

        /* Check detection of argument 3   */
        M = -1; GESTEIN_CALL(3, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M = 0;
        /* Check detection of argument 5   */
        LDA = 0; GESTEIN_CALL(5, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M = 2; LDA=1; GESTEIN_CALL(5, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDA = 5; M = 0;

        /* Check detection of argument 7   */
        LDQ = 0; GESTEIN_CALL(7, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M = 2; LDQ=1; GESTEIN_CALL(7, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDQ = 5; M = 0;

        /* Check detection of argument 9   */
        LDX = 0; GESTEIN_CALL(9, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M = 2; LDX=1; GESTEIN_CALL(9, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDX = 5; M = 0;

        /* Check detection of argument 12   */
        LDWORK = 0; GESTEIN_CALL(12, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDWORK = 5000000;


        /* Check Rest  */
        if ( trstein_functions[i].chk_l3 ) {
            mepack_single_trstein_blocksize_set(1);
            GESTEIN_CALL(30, "F" ,"N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_single_trstein_blocksize_set(2);
        }
        if ( trstein_functions[i].isolver == MEPACK_FRONTEND_SOLVER_2STAGE ) {
            mepack_single_trstein_blocksize_set(4);
            mepack_single_trstein_blocksize_2stage_set(2);
            GESTEIN_CALL(33, "F", "N", &M, NULL, &LDA, NULL, &LDQ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_single_trstein_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

