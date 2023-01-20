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
typedef void (*trlyap_function) (const char*FACT, const char *TRANSA, Int *M, double * A, Int *LDA, double *B, Int *LDB,
        double *Q, Int *LDQ, double *Z, Int * LDZ, double *X, Int *LDX, double * SCALE, double *WORK, Int *LDWORK,
        Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


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

#define GELYAP_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17)  do { cname = trlyap_functions[i].name; global_state =0;  \
            ign_wrng = trlyap_functions[i].ign_wrng_name, count_all++;  \
            trlyap(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,1, 1); \
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
    double WORK = 0 , SCALE = 0;
    Int M = 0, LDA=0, LDX=0, INFO=0, LDQ = 0;
    Int LDZ = 0, LDB = 0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    Int LDWORK;
    trlyap_function trlyap;
    struct function_name {
        trlyap_function fn;
        char * name;
        int chk_l3;
        int isolver;
        int ign_wrng_name;
    };

    struct function_name trlyap_functions[] = {
        { &(FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)), "DLA_GGSTEIN", 1, MEPACK_FRONTEND_SOLVER_LEVEL3 ,0 },
        { &(FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)), "DLA_GGSTEIN", 1, MEPACK_FRONTEND_SOLVER_LEVEL2 ,0 },
        { &(FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)), "DLA_GGSTEIN", 1, MEPACK_FRONTEND_SOLVER_DAG ,0 },
        { &(FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)), "DLA_GGSTEIN", 1, MEPACK_FRONTEND_SOLVER_2STAGE ,0 },
        { &(FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)), "DLA_GGSTEIN", 1, MEPACK_FRONTEND_SOLVER_RECURSIVE ,0 },

        { NULL, NULL, 0,0, 0}
    };

    argc = argc;
    (void) argv;

    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( trlyap_functions[i].fn != NULL ) {
        printf("Check %s\n", trlyap_functions[i].name);
        LDWORK = 1024*1024;
        mepack_tgstein_frontend_solver_set(trlyap_functions[i].isolver);

        trlyap = trlyap_functions[i].fn;
         /* Check detection of argument 1   */
        GELYAP_CALL(1, "B", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        /* Check detection of argument 2   */
        GELYAP_CALL(2, "N", "X", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        /* Check detection of argument 3   */
        M=-1; GELYAP_CALL(3, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M = 0;

        /* Check detection of argument 5   */
        LDA=0; GELYAP_CALL(5, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M=2; LDA=1; GELYAP_CALL(5, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDA = 5;

        /* Check detection of argument 7   */
        LDB=0; GELYAP_CALL(7, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M=2; LDB=1; GELYAP_CALL(7, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDB = 5;

        /* Check detection of argument 9   */
        LDQ=0; GELYAP_CALL(9, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M=2; LDQ=1; GELYAP_CALL(9, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDQ = 5;

        /* Check detection of argument 11   */
        LDZ=0; GELYAP_CALL(11, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M=2; LDZ=1; GELYAP_CALL(11, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDZ = 5;

        /* Check detection of argument 13   */
        LDX=0; GELYAP_CALL(13, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        M=2; LDX=1; GELYAP_CALL(13, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDX = 5;

        /* Check detection of argument 16   */
        LDWORK = 0; GELYAP_CALL(16, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
        LDWORK = 1024*1024;

        if ( trlyap_functions[i].chk_l3 ) {
            mepack_double_tgstein_blocksize_set(1);
            GELYAP_CALL(30, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_double_tgstein_blocksize_set(2);

        }
        if ( trlyap_functions[i].isolver == MEPACK_FRONTEND_SOLVER_2STAGE ) {
            mepack_double_tgstein_blocksize_set(4);
            mepack_double_tgstein_blocksize_2stage_set(2);
            GELYAP_CALL(33, "N", "N", &M, NULL, &LDA, NULL,  &LDB, NULL, &LDQ, NULL, &LDZ, NULL, &LDX, &SCALE, &WORK, &LDWORK, &INFO);
            mepack_double_tgstein_blocksize_2stage_set(512);
        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

