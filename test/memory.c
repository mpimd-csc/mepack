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
#include <unistd.h>
#include <stdint.h>


#include "mepack.h"

struct memory_functions_st {
    char* c_function;
    char* f_function;
} memfunctions[] = {
   {"mepack_double_tgcsylv_dag", "DLA_TGCSYLV_DAG"},
   {"mepack_double_tgcsylv_dual_dag", "DLA_TGCSYLV_DUAL_DAG"},
   {"mepack_double_tgcsylv_dual_level2", "DLA_TGCSYLV_DUAL_L2"},
   {"mepack_double_tgcsylv_dual_level2_local_copy", "DLA_TGCSYLV_DUAL_L2_LOCAL_COPY"},
   {"mepack_double_tgcsylv_dual_level2_local_copy_128", "DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128"},
   {"mepack_double_tgcsylv_dual_level2_local_copy_32", "DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32"},
   {"mepack_double_tgcsylv_dual_level2_local_copy_64", "DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64"},
   {"mepack_double_tgcsylv_dual_level2_local_copy_96", "DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96"},
   {"mepack_double_tgcsylv_dual_level3", "DLA_TGCSYLV_DUAL_L3"},
   {"mepack_double_tgcsylv_dual_level3_2stage", "DLA_TGCSYLV_DUAL_L3_2S"},
   {"mepack_double_tgcsylv_dual_recursive", "DLA_TGCSYLV_DUAL_RECURSIVE"},
   {"mepack_double_tgcsylv_level2", "DLA_TGCSYLV_L2"},
   {"mepack_double_tgcsylv_level2_local_copy", "DLA_TGCSYLV_L2_LOCAL_COPY"},
   {"mepack_double_tgcsylv_level2_local_copy_128", "DLA_TGCSYLV_L2_LOCAL_COPY_128"},
   {"mepack_double_tgcsylv_level2_local_copy_32", "DLA_TGCSYLV_L2_LOCAL_COPY_32"},
   {"mepack_double_tgcsylv_level2_local_copy_64", "DLA_TGCSYLV_L2_LOCAL_COPY_64"},
   {"mepack_double_tgcsylv_level2_local_copy_96", "DLA_TGCSYLV_L2_LOCAL_COPY_96"},
   {"mepack_double_tgcsylv_level2_reorder", "DLA_TGCSYLV_L2_REORDER"},
   {"mepack_double_tgcsylv_level2_unopt", "DLA_TGCSYLV_L2_UNOPT"},
   {"mepack_double_tgcsylv_level3", "DLA_TGCSYLV_L3"},
   {"mepack_double_tgcsylv_level3_2stage", "DLA_TGCSYLV_L3_2S"},
   {"mepack_double_tgcsylv_level3_unopt", "DLA_TGCSYLV_L3_UNOPT"},
   {"mepack_double_tgcsylv_recursive", "DLA_TGCSYLV_RECURSIVE"},
   {"mepack_double_tglyap_dag", "DLA_TGLYAP_DAG"},
   {"mepack_double_tglyap_level2", "DLA_TGLYAP_L2"},
   {"mepack_double_tglyap_level2_unopt", "DLA_TGLYAP_L2_UNOPT"},
   {"mepack_double_tglyap_level3", "DLA_TGLYAP_L3"},
   {"mepack_double_tglyap_level3_2stage", "DLA_TGLYAP_L3_2S"},
   {"mepack_double_tglyap_recursive", "DLA_TGLYAP_RECURSIVE"},
   {"mepack_double_tgstein_dag", "DLA_TGSTEIN_DAG"},
   {"mepack_double_tgstein_level2", "DLA_TGSTEIN_L2"},
   {"mepack_double_tgstein_level3", "DLA_TGSTEIN_L3"},
   {"mepack_double_tgstein_level3_2stage", "DLA_TGSTEIN_L3_2S"},
   {"mepack_double_tgstein_recursive", "DLA_TGSTEIN_RECURSIVE"},
   {"mepack_double_tgsylv_dag", "DLA_TGSYLV_DAG"},
   {"mepack_double_tgsylv_gardiner_laub", "DLA_TGSYLV_GARDINER_LAUB"},
   {"mepack_double_tgsylv_level2", "DLA_TGSYLV_L2"},
   {"mepack_double_tgsylv_level2_colwise", "DLA_TGSYLV_L2_COLWISE"},
   {"mepack_double_tgsylv_level2_local_copy", "DLA_TGSYLV_L2_LOCAL_COPY"},
   {"mepack_double_tgsylv_level2_local_copy_128", "DLA_TGSYLV_L2_LOCAL_COPY_128"},
   {"mepack_double_tgsylv_level2_local_copy_32", "DLA_TGSYLV_L2_LOCAL_COPY_32"},
   {"mepack_double_tgsylv_level2_local_copy_64", "DLA_TGSYLV_L2_LOCAL_COPY_64"},
   {"mepack_double_tgsylv_level2_local_copy_96", "DLA_TGSYLV_L2_LOCAL_COPY_96"},
   {"mepack_double_tgsylv_level2_reorder", "DLA_TGSYLV_L2_REORDER"},
   {"mepack_double_tgsylv_level2_rowwise", "DLA_TGSYLV_L2_ROWWISE"},
   {"mepack_double_tgsylv_level3", "DLA_TGSYLV_L3"},
   {"mepack_double_tgsylv_level3_2stage", "DLA_TGSYLV_L3_2S"},
   {"mepack_double_tgsylv_level3_colwise", "DLA_TGSYLV_L3_COLWISE"},
   {"mepack_double_tgsylv_level3_rowwise", "DLA_TGSYLV_L3_ROWWISE"},
   {"mepack_double_tgsylv_recursive", "DLA_TGSYLV_RECURSIVE"},
   {"mepack_double_trlyap_dag", "DLA_TRLYAP_DAG"},
   {"mepack_double_trlyap_level2", "DLA_TRLYAP_L2"},
   {"mepack_double_trlyap_level2_opt", "DLA_TRLYAP_L2_OPT"},
   {"mepack_double_trlyap_level3", "DLA_TRLYAP_L3"},
   {"mepack_double_trlyap_level3_2stage", "DLA_TRLYAP_L3_2S"},
   {"mepack_double_trlyap_recursive", "DLA_TRLYAP_RECURSIVE"},
   {"mepack_double_trstein_dag", "DLA_TRSTEIN_DAG"},
   {"mepack_double_trstein_level2", "DLA_TRSTEIN_L2"},
   {"mepack_double_trstein_level3", "DLA_TRSTEIN_L3"},
   {"mepack_double_trstein_level3_2stage", "DLA_TRSTEIN_L3_2S"},
   {"mepack_double_trstein_recursive", "DLA_TRSTEIN_RECURSIVE"},
   {"mepack_double_trsylv2_dag", "DLA_TRSYLV2_DAG"},
   {"mepack_double_trsylv2_level2", "DLA_TRSYLV2_L2"},
   {"mepack_double_trsylv2_level2_local_copy", "DLA_TRSYLV2_L2_LOCAL_COPY"},
   {"mepack_double_trsylv2_level2_local_copy_128", "DLA_TRSYLV2_L2_LOCAL_COPY_128"},
   {"mepack_double_trsylv2_level2_local_copy_32", "DLA_TRSYLV2_L2_LOCAL_COPY_32"},
   {"mepack_double_trsylv2_level2_local_copy_64", "DLA_TRSYLV2_L2_LOCAL_COPY_64"},
   {"mepack_double_trsylv2_level2_local_copy_96", "DLA_TRSYLV2_L2_LOCAL_COPY_96"},
   {"mepack_double_trsylv2_level2_reorder", "DLA_TRSYLV2_L2_REORDER"},
   {"mepack_double_trsylv2_level2_unopt", "DLA_TRSYLV2_L2_UNOPT"},
   {"mepack_double_trsylv2_level3", "DLA_TRSYLV2_L3"},
   {"mepack_double_trsylv2_level3_2stage", "DLA_TRSYLV2_L3_2S"},
   {"mepack_double_trsylv2_level3_unopt", "DLA_TRSYLV2_L3_UNOPT"},
   {"mepack_double_trsylv2_recursive", "DLA_TRSYLV2_RECURSIVE"},
   {"mepack_double_trsylv_dag", "DLA_TRSYLV_DAG"},
   {"mepack_double_trsylv_level2", "DLA_TRSYLV_L2"},
   {"mepack_double_trsylv_level2_local_copy", "DLA_TRSYLV_L2_LOCAL_COPY"},
   {"mepack_double_trsylv_level2_local_copy_128", "DLA_TRSYLV_L2_LOCAL_COPY_128"},
   {"mepack_double_trsylv_level2_local_copy_32", "DLA_TRSYLV_L2_LOCAL_COPY_32"},
   {"mepack_double_trsylv_level2_local_copy_64", "DLA_TRSYLV_L2_LOCAL_COPY_64"},
   {"mepack_double_trsylv_level2_local_copy_96", "DLA_TRSYLV_L2_LOCAL_COPY_96"},
   {"mepack_double_trsylv_level2_reorder", "DLA_TRSYLV_L2_REORDER"},
   {"mepack_double_trsylv_level2_unopt", "DLA_TRSYLV_L2_UNOPT"},
   {"mepack_double_trsylv_level3", "DLA_TRSYLV_L3"},
   {"mepack_double_trsylv_level3_2stage", "DLA_TRSYLV_L3_2S"},
   {"mepack_double_trsylv_level3_unopt", "DLA_TRSYLV_L3_UNOPT"},
   {"mepack_double_trsylv_recursive", "DLA_TRSYLV_RECURSIVE"},
   {"mepack_single_tgcsylv_dag", "SLA_TGCSYLV_DAG"},
   {"mepack_single_tgcsylv_dual_dag", "SLA_TGCSYLV_DUAL_DAG"},
   {"mepack_single_tgcsylv_dual_level2", "SLA_TGCSYLV_DUAL_L2"},
   {"mepack_single_tgcsylv_dual_level2_local_copy", "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY"},
   {"mepack_single_tgcsylv_dual_level2_local_copy_128", "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128"},
   {"mepack_single_tgcsylv_dual_level2_local_copy_32", "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32"},
   {"mepack_single_tgcsylv_dual_level2_local_copy_64", "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64"},
   {"mepack_single_tgcsylv_dual_level2_local_copy_96", "SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96"},
   {"mepack_single_tgcsylv_dual_level3", "SLA_TGCSYLV_DUAL_L3"},
   {"mepack_single_tgcsylv_dual_level3_2stage", "SLA_TGCSYLV_DUAL_L3_2S"},
   {"mepack_single_tgcsylv_dual_recursive", "SLA_TGCSYLV_DUAL_RECURSIVE"},
   {"mepack_single_tgcsylv_level2", "SLA_TGCSYLV_L2"},
   {"mepack_single_tgcsylv_level2_local_copy", "SLA_TGCSYLV_L2_LOCAL_COPY"},
   {"mepack_single_tgcsylv_level2_local_copy_128", "SLA_TGCSYLV_L2_LOCAL_COPY_128"},
   {"mepack_single_tgcsylv_level2_local_copy_32", "SLA_TGCSYLV_L2_LOCAL_COPY_32"},
   {"mepack_single_tgcsylv_level2_local_copy_64", "SLA_TGCSYLV_L2_LOCAL_COPY_64"},
   {"mepack_single_tgcsylv_level2_local_copy_96", "SLA_TGCSYLV_L2_LOCAL_COPY_96"},
   {"mepack_single_tgcsylv_level2_reorder", "SLA_TGCSYLV_L2_REORDER"},
   {"mepack_single_tgcsylv_level2_unopt", "SLA_TGCSYLV_L2_UNOPT"},
   {"mepack_single_tgcsylv_level3", "SLA_TGCSYLV_L3"},
   {"mepack_single_tgcsylv_level3_2stage", "SLA_TGCSYLV_L3_2S"},
   {"mepack_single_tgcsylv_level3_unopt", "SLA_TGCSYLV_L3_UNOPT"},
   {"mepack_single_tgcsylv_recursive", "SLA_TGCSYLV_RECURSIVE"},
   {"mepack_single_tglyap_dag", "SLA_TGLYAP_DAG"},
   {"mepack_single_tglyap_level2", "SLA_TGLYAP_L2"},
   {"mepack_single_tglyap_level2_unopt", "SLA_TGLYAP_L2_UNOPT"},
   {"mepack_single_tglyap_level3", "SLA_TGLYAP_L3"},
   {"mepack_single_tglyap_level3_2stage", "SLA_TGLYAP_L3_2S"},
   {"mepack_single_tglyap_recursive", "SLA_TGLYAP_RECURSIVE"},
   {"mepack_single_tgstein_dag", "SLA_TGSTEIN_DAG"},
   {"mepack_single_tgstein_level2", "SLA_TGSTEIN_L2"},
   {"mepack_single_tgstein_level3", "SLA_TGSTEIN_L3"},
   {"mepack_single_tgstein_level3_2stage", "SLA_TGSTEIN_L3_2S"},
   {"mepack_single_tgstein_recursive", "SLA_TGSTEIN_RECURSIVE"},
   {"mepack_single_tgsylv_dag", "SLA_TGSYLV_DAG"},
   {"mepack_single_tgsylv_gardiner_laub", "SLA_TGSYLV_GARDINER_LAUB"},
   {"mepack_single_tgsylv_level2", "SLA_TGSYLV_L2"},
   {"mepack_single_tgsylv_level2_colwise", "SLA_TGSYLV_L2_COLWISE"},
   {"mepack_single_tgsylv_level2_local_copy", "SLA_TGSYLV_L2_LOCAL_COPY"},
   {"mepack_single_tgsylv_level2_local_copy_128", "SLA_TGSYLV_L2_LOCAL_COPY_128"},
   {"mepack_single_tgsylv_level2_local_copy_32", "SLA_TGSYLV_L2_LOCAL_COPY_32"},
   {"mepack_single_tgsylv_level2_local_copy_64", "SLA_TGSYLV_L2_LOCAL_COPY_64"},
   {"mepack_single_tgsylv_level2_local_copy_96", "SLA_TGSYLV_L2_LOCAL_COPY_96"},
   {"mepack_single_tgsylv_level2_reorder", "SLA_TGSYLV_L2_REORDER"},
   {"mepack_single_tgsylv_level2_rowwise", "SLA_TGSYLV_L2_ROWWISE"},
   {"mepack_single_tgsylv_level3", "SLA_TGSYLV_L3"},
   {"mepack_single_tgsylv_level3_2stage", "SLA_TGSYLV_L3_2S"},
   {"mepack_single_tgsylv_level3_colwise", "SLA_TGSYLV_L3_COLWISE"},
   {"mepack_single_tgsylv_level3_rowwise", "SLA_TGSYLV_L3_ROWWISE"},
   {"mepack_single_tgsylv_recursive", "SLA_TGSYLV_RECURSIVE"},
   {"mepack_single_trlyap_dag", "SLA_TRLYAP_DAG"},
   {"mepack_single_trlyap_level2", "SLA_TRLYAP_L2"},
   {"mepack_single_trlyap_level2_opt", "SLA_TRLYAP_L2_OPT"},
   {"mepack_single_trlyap_level3", "SLA_TRLYAP_L3"},
   {"mepack_single_trlyap_level3_2stage", "SLA_TRLYAP_L3_2S"},
   {"mepack_single_trlyap_recursive", "SLA_TRLYAP_RECURSIVE"},
   {"mepack_single_trstein_dag", "SLA_TRSTEIN_DAG"},
   {"mepack_single_trstein_level2", "SLA_TRSTEIN_L2"},
   {"mepack_single_trstein_level3", "SLA_TRSTEIN_L3"},
   {"mepack_single_trstein_level3_2stage", "SLA_TRSTEIN_L3_2S"},
   {"mepack_single_trstein_recursive", "SLA_TRSTEIN_RECURSIVE"},
   {"mepack_single_trsylv2_dag", "SLA_TRSYLV2_DAG"},
   {"mepack_single_trsylv2_level2", "SLA_TRSYLV2_L2"},
   {"mepack_single_trsylv2_level2_local_copy", "SLA_TRSYLV2_L2_LOCAL_COPY"},
   {"mepack_single_trsylv2_level2_local_copy_128", "SLA_TRSYLV2_L2_LOCAL_COPY_128"},
   {"mepack_single_trsylv2_level2_local_copy_32", "SLA_TRSYLV2_L2_LOCAL_COPY_32"},
   {"mepack_single_trsylv2_level2_local_copy_64", "SLA_TRSYLV2_L2_LOCAL_COPY_64"},
   {"mepack_single_trsylv2_level2_local_copy_96", "SLA_TRSYLV2_L2_LOCAL_COPY_96"},
   {"mepack_single_trsylv2_level2_reorder", "SLA_TRSYLV2_L2_REORDER"},
   {"mepack_single_trsylv2_level2_unopt", "SLA_TRSYLV2_L2_UNOPT"},
   {"mepack_single_trsylv2_level3", "SLA_TRSYLV2_L3"},
   {"mepack_single_trsylv2_level3_2stage", "SLA_TRSYLV2_L3_2S"},
   {"mepack_single_trsylv2_level3_unopt", "SLA_TRSYLV2_L3_UNOPT"},
   {"mepack_single_trsylv2_recursive", "SLA_TRSYLV2_RECURSIVE"},
   {"mepack_single_trsylv_dag", "SLA_TRSYLV_DAG"},
   {"mepack_single_trsylv_level2", "SLA_TRSYLV_L2"},
   {"mepack_single_trsylv_level2_local_copy", "SLA_TRSYLV_L2_LOCAL_COPY"},
   {"mepack_single_trsylv_level2_local_copy_128", "SLA_TRSYLV_L2_LOCAL_COPY_128"},
   {"mepack_single_trsylv_level2_local_copy_32", "SLA_TRSYLV_L2_LOCAL_COPY_32"},
   {"mepack_single_trsylv_level2_local_copy_64", "SLA_TRSYLV_L2_LOCAL_COPY_64"},
   {"mepack_single_trsylv_level2_local_copy_96", "SLA_TRSYLV_L2_LOCAL_COPY_96"},
   {"mepack_single_trsylv_level2_reorder", "SLA_TRSYLV_L2_REORDER"},
   {"mepack_single_trsylv_level2_unopt", "SLA_TRSYLV_L2_UNOPT"},
   {"mepack_single_trsylv_level3", "SLA_TRSYLV_L3"},
   {"mepack_single_trsylv_level3_2stage", "SLA_TRSYLV_L3_2S"},
   {"mepack_single_trsylv_level3_unopt", "SLA_TRSYLV_L3_UNOPT"},
   {"mepack_single_trsylv_recursive", "SLA_TRSYLV_RECURSIVE"},
   {NULL, NULL}
};

struct memory_functions_st memfunctions_frontend[] =
{
    {"mepack_double_gelyap", "DLA_GELYAP"},
    {"mepack_single_gelyap", "SLA_GELYAP"},
    {"mepack_double_gestein", "DLA_GESTEIN"},
    {"mepack_single_gestein", "SLA_GESTEIN"},
    {"mepack_double_gesylv", "DLA_GESYLV"},
    {"mepack_single_gesylv", "SLA_GESYLV"},
    {"mepack_double_gesylv2", "DLA_GESYLV2"},
    {"mepack_single_gesylv2", "SLA_GESYLV2"},
    {"mepack_double_gglyap", "DLA_GGLYAP"},
    {"mepack_single_gglyap", "SLA_GGLYAP"},
    {"mepack_double_ggstein", "DLA_GGSTEIN"},
    {"mepack_single_ggstein", "SLA_GGSTEIN"},
    {"mepack_double_ggsylv", "DLA_GGSYLV"},
    {"mepack_single_ggsylv", "SLA_GGSYLV"},
    {"mepack_double_ggcsylv", "DLA_GGCSYLV"},
    {"mepack_single_ggcsylv", "SLA_GGCSYLV"},
    {"mepack_double_ggcsylv_dual", "DLA_GGCSYLV_DUAL"},
    {"mepack_single_ggcsylv_dual", "SLA_GGCSYLV_DUAL"},
    {"mepack_double_gelyap_refine", "DLA_GELYAP_REFINE"},
    {"mepack_single_gelyap_refine", "SLA_GELYAP_REFINE"},
    {"mepack_double_gestein_refine", "DLA_GESTEIN_REFINE"},
    {"mepack_single_gestein_refine", "SLA_GESTEIN_REFINE"},
    {"mepack_double_gglyap_refine", "DLA_GGLYAP_REFINE"},
    {"mepack_single_gglyap_refine", "SLA_GGLYAP_REFINE"},
    {"mepack_double_ggstein_refine", "DLA_GGSTEIN_REFINE"},
    {"mepack_single_ggstein_refine", "SLA_GGSTEIN_REFINE"},
    {"mepack_double_gesylv_refine", "DLA_GESYLV_REFINE"},
    {"mepack_single_gesylv_refine", "SLA_GESYLV_REFINE"},
    {"mepack_double_gesylv2_refine", "DLA_GESYLV2_REFINE"},
    {"mepack_single_gesylv2_refine", "SLA_GESYLV2_REFINE"},
    {"mepack_double_ggsylv_refine", "DLA_GGSYLV_REFINE"},
    {"mepack_single_ggsylv_refine", "SLA_GGSYLV_REFINE"},
    {"mepack_double_ggcsylv_refine", "DLA_GGCSYLV_REFINE"},
    {"mepack_single_ggcsylv_refine", "SLA_GGCSYLV_REFINE"},
    {"mepack_double_ggcsylv_dual_refine", "DLA_GGCSYLV_DUAL_REFINE"},
    {"mepack_single_ggcsylv_dual_refine", "SLA_GGCSYLV_DUAL_REFINE"},
    {NULL, NULL}
};



int Ms[] = { 0, 1,  94, 257, 519, 1901, -1};
int Ns[] = { 0, 1, 3, 65, 367, 753, -1};

struct fact_st {
    char *factA;
    char *factB;
} facts [] = {
    {"N", "N"},
    {"N", "F"},
    {"N", "H"},
    {"F", "N"},
    {"F", "F"},
    {"F", "H"},
    {"H", "N"},
    {"H", "F"},
    {"H", "H"},
    {NULL, NULL}
};


int test_memory(void)
{
    int count = 0, fn, k, l;
    ssize_t r1 , r2;

    fn = 0;
    while ( memfunctions[fn].c_function != NULL)
    {
        k = 0;
        while (Ms[k] >= 0)
        {
            l = 0;
            while(Ns[l] >= 0)
            {
                r1 = mepack_memory(memfunctions[fn].c_function, Ms[k], Ns[l]);
                r2 = mepack_memory(memfunctions[fn].f_function, Ms[k], Ns[l]);

                fprintf(stdout, "mem(%s,%s) (%d, %d)  = (%ld, %ld) %s\n", memfunctions[fn].c_function, memfunctions[fn].f_function, Ms[k], Ns[l],
                        r1, r2, (r1 == r2) ? "OK" : "Error");
                if ( r1 != r2)  count++;
                l++;
            }
            k++;
        }
        fn++;
    }

    return count;

}


int test_memory_frontend(void)
{
    int count = 0, fn, k, l, j ;
    ssize_t r1 , r2;

    fn = 0;
    while ( memfunctions_frontend[fn].c_function != NULL)
    {
        k = 0;
        while (Ms[k] >= 0)
        {
            l = 0;
            while(Ns[l] >= 0)
            {
                j = 0;
                while( facts[j].factA != NULL)
                {
                    r1 = mepack_memory_frontend(memfunctions_frontend[fn].c_function, facts[j].factA, facts[j].factB, Ms[k], Ns[l]);
                    r2 = mepack_memory_frontend(memfunctions_frontend[fn].f_function, facts[j].factA, facts[j].factB, Ms[k], Ns[l]);

                    fprintf(stdout, "mem(%s,%s) (%s,%s) (%d, %d)  = (%ld, %ld) %s\n", memfunctions_frontend[fn].c_function, memfunctions_frontend[fn].f_function,
                            facts[j].factA, facts[j].factB, Ms[k], Ns[l],
                            r1, r2, (r1 == r2) ? "OK" : "Error");
                    if ( r1 != r2)  count++;
                    j++;
                }
                l++;
            }
            k++;
        }
        fn++;
    }

    return count;


    return 0;
}



int main(int argc, char ** argv)
{
    (void) argc;
    (void) argv;
    int ret1 = 0;
    int ret2 = 0;

    /** ret1 = test_memory(); */

    ret2 = test_memory_frontend();

    return ret1 + ret2;

}

