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
#include <string.h>
#include <complex.h>
#include <math.h>
#include "proto.h"

typedef void (*tgcsylv_function) (const char *TRANSA, const char*TRANSB, float *sgn, float * sgn2,  Int *M, Int *N, float * A, Int *LDA,float * B, Int * LDB, float *C, Int *ldc,float* D, Int *LDC, float *X, Int *LDX, float *Y, Int * LDY,
        float * SCALE, float *WORK, Int *INFO, FORTRAN_LEN_T l1, FORTRAN_LEN_T l2);


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

#define TGCSYLV_CALL(D, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21)  do { cname = tgcsylv_functions[i].name; global_state =0 ;  \
            ign_wrng = tgcsylv_functions[i].ign_wrng_name, count_all++;  \
            tgcsylv(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,1,1); \
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
    Int LDY = 0;
    float sign2 = 0;
    int ok = 1;
    int i = 0;
    int count_ok = 0;
    int count_failed = 0;
    int count_all = 0;
    tgcsylv_function tgcsylv;
    struct function_name {
        tgcsylv_function fn;
        char * name;
        int chk_l3;
        int chk_l3_2s;
        int ign_wrng_name;
        int dual;
    };

    struct function_name tgcsylv_functions[] = {
        { &(FC_GLOBAL_(sla_tgcsylv_l3,SLA_TGCSYLV_L3)), "SLA_TGCSYLV_L3", 1, 0, 0, 0 },
        { &(FC_GLOBAL_(sla_tgcsylv_l3_unopt,SLA_TGCSYLV_L3_UNOPT)), "SLA_TGCSYLV_L3_UNOPT", 1, 0,0,0  },
        { &(FC_GLOBAL_(sla_tgcsylv_l2,SLA_TGCSYLV_L2)), "SLA_TGCSYLV_L2", 0, 0 ,0,0 },
        { &(FC_GLOBAL_(sla_tgcsylv_l2_reorder,SLA_TGCSYLV_L2_REORDER)), "SLA_TGCSYLV_L2_REORDER", 0, 0,0,0  },
        { &(FC_GLOBAL_(sla_tgcsylv_l2_local_copy,SLA_TGCSYLV_L2_LOCAL_COPY)), "SLA_TGCSYLV_L2_LOCAL_COPY" ,0 ,0,0,0  },
        { &(FC_GLOBAL_(sla_tgcsylv_l2_local_copy_32,SLA_TGCSYLV_L2_LOCAL_COPY_32)), "SLA_TGCSYLV_L2_LOCAL_COPY_32", 0,0 , 1,0},
        { &(FC_GLOBAL_(sla_tgcsylv_l2_local_copy_64,SLA_TGCSYLV_L2_LOCAL_COPY_64)), "SLA_TGCSYLV_L2_LOCAL_COPY_64", 0,0 , 1,0},
        { &(FC_GLOBAL_(sla_tgcsylv_l2_local_copy_96,SLA_TGCSYLV_L2_LOCAL_COPY_96)), "SLA_TGCSYLV_L2_LOCAL_COPY_96", 0,0 , 1,0},
        { &(FC_GLOBAL_(sla_tgcsylv_l2_local_copy_128,SLA_TGCSYLV_L2_LOCAL_COPY_128)), "SLA_TGCSYLV_L2_LOCAL_COPY_128", 0,0,1,0 },
        { &(FC_GLOBAL_(sla_tgcsylv_l2_unopt,SLA_TGCSYLV_L2_UNOPT)), "SLA_TGCSYLV_L2_UNOPT", 0,0, 0,0 },
        { &(FC_GLOBAL_(sla_tgcsylv_dag,SLA_TGCSYLV_DAG)), "SLA_TGCSYLV_DAG", 1, 0, 0,0 },
        { &(FC_GLOBAL_(sla_tgcsylv_l3_2s,SLA_TGCSYLV_L3_2S)), "SLA_TGCSYLV_L3_2S",1 , 1 ,0,0},
        { &(FC_GLOBAL_(sla_tgcsylv_recursive,SLA_TGCSYLV_RECURSIVE)), "SLA_TGCSYLV_RECURSIVE", 0, 0 , 0,0},
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l3,SLA_TGCSYLV_DUAL_L3)), "SLA_TGCSYLV_DUAL_L3", 1, 0,0,1 },
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l2,SLA_TGCSYLV_DUAL_L2)), "SLA_TGCSYLV_DUAL_L2", 0, 0 ,0,1 },
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY)), "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY" ,0 ,0,0,1 },
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_32,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)), "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32", 0,0 , 1,1},
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_64,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)), "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64", 0,0 , 1,1},
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_96,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)), "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96", 0,0 , 1,1},
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_128,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)), "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128", 0,0,1,1 },
        { &(FC_GLOBAL_(sla_tgcsylv_dual_dag,SLA_TGCSYLV_DUAL_DAG)), "SLA_TGCSYLV_DUAL_DAG", 1, 0, 0,1 },
        { &(FC_GLOBAL_(sla_tgcsylv_dual_l3_2s,SLA_TGCSYLV_DUAL_L3_2S)), "SLA_TGCSYLV_DUAL_L3_2S",1 , 1 ,0,1},
        { &(FC_GLOBAL_(sla_tgcsylv_dual_recursive,SLA_TGCSYLV_DUAL_RECURSIVE)), "SLA_TGCSYLV_DUAL_RECURSIVE", 0, 0 , 0,1},

        { NULL, NULL, 0,0, 0, 0}
    };


    argc = argc;
    (void) argv;


    /* Initialize parameters  */
    /* Initialize Error Handler  */
    FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C)(err_handler);

    while ( tgcsylv_functions[i].fn != NULL ) {
        printf("Check %s\n", tgcsylv_functions[i].name);
        tgcsylv = tgcsylv_functions[i].fn;

        /* Check detection of argument 1   */
        TGCSYLV_CALL(1, "X", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        /* Check detection of argument 2   */
        TGCSYLV_CALL(2, "N", "X", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        /* Check detection of argument 3   */
        sign = 0;   TGCSYLV_CALL(3, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO );
        sign = 2;   TGCSYLV_CALL(3, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        sign = -2;  TGCSYLV_CALL(3, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        sign = 1;
         /* Check detection of argument 4   */
        sign2 = 0;   TGCSYLV_CALL(4, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        sign2 = 2;   TGCSYLV_CALL(4, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        sign2 = -2;  TGCSYLV_CALL(4, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        sign2 = 1;

        /* Check detection of argument 5  */
        M = -1; TGCSYLV_CALL(5, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 0;
        /* Check detection of argument 6  */
        N = -1; TGCSYLV_CALL(6, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        N = 0;
        /* Check detection of argument 8  */
        LDA = 0;        TGCSYLV_CALL(8, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 2; LDA = 1; TGCSYLV_CALL(8, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 0; LDA = 5;
        /* Check detection of argument 10  */
        LDB = 0; TGCSYLV_CALL(10, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        N = 2; LDB = 1; TGCSYLV_CALL(10, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        N = 0; LDB = 5;

        /* Check detection of argument 12  */
        LDC = 0;        TGCSYLV_CALL(12, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 2; LDC = 1; TGCSYLV_CALL(12, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 0; LDC = 5;

        /* Check detection of argument 14  */
        LDD = 0; TGCSYLV_CALL(14, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        N = 2; LDD = 1; TGCSYLV_CALL(14, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        N = 0; LDD = 5;

        /* Check detection of argument 16  */
        LDX = 0; TGCSYLV_CALL(16, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 2; LDA = 2;  LDX = 1; TGCSYLV_CALL(16, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 0; LDA = 1;  LDX = 5;

        /* Check detection of argument 18  */
        LDY = 0; TGCSYLV_CALL(18, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB, NULL, &LDC, NULL, &LDD, NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 2; LDA = 2;  LDY = 1; TGCSYLV_CALL(18, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
        M = 0; LDA = 1;  LDY = 5;


        if ( tgcsylv_functions[i].chk_l3 ) {
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_mb_set(1); else mepack_single_tgcsylv_blocksize_mb_set(1);
            TGCSYLV_CALL(30, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO );
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_mb_set(2); else mepack_single_tgcsylv_blocksize_mb_set(2);

            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_nb_set(1); else mepack_single_tgcsylv_blocksize_nb_set(1);
            TGCSYLV_CALL(31, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_nb_set(2); else mepack_single_tgcsylv_blocksize_nb_set(2);
        }
        if ( tgcsylv_functions[i].chk_l3_2s ) {
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_mb_set(4); else mepack_single_tgcsylv_blocksize_mb_set(4);
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_nb_set(4); else mepack_single_tgcsylv_blocksize_nb_set(4);
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_2stage_set(2); else mepack_single_tgcsylv_blocksize_2stage_set(2);

            TGCSYLV_CALL(33, "N", "N", &sign, &sign2, &M, &N, NULL, &LDA, NULL, &LDB,NULL, &LDC, NULL, &LDD,  NULL, &LDX, NULL, &LDY, &SCALE, &WORK, &INFO);
            if (tgcsylv_functions[i].dual) mepack_single_tgcsylv_dual_blocksize_2stage_set(512); else mepack_single_tgcsylv_blocksize_2stage_set(512);

        }

        printf("\n");
        i++;
    }

    printf("Tests: (%d/%d/%d)\n", count_all, count_ok, count_failed);
    return (count_all==count_ok)? 0:-1;
}

