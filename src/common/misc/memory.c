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


#include "FCMangle.h"
#include "mepack.h"
#include "mepack_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX(A,B) ((A)>(B)?(A):(B))

/**
 * @addtogroup auxmem
 * @{
 */

/**
 * \brief Query the size of the workspace for triangular solvers.
 *
 * \par Purpose:
 *
 * The mepack_memory function queries the size of the required routine. It wraps around the Fortran
 * style workspace query and adds a typecast to avoid integer overflows. This can happen if MEPACK
 * is compiled with the support for 64-bit integers in Fortran. If solver only relies on one dimension
 * argument, e.g. a Lyapunov or a Stein equation is planed to be solved, only the \b M parameter is used.
 * The value of \b N is neglected in this case. The frontend routines, i.e. the ones dealing with the
 * transformation matrices as well, have to use \ref mepack_memory_frontend to query their memory
 * requirements.
 *
 * \see mepack_memory_frontend
 *
 * \remark
 *    This function should be used from C. If the workspace query is performed from a Fortran code,
 *    the workspace is queried by passing `INFO == -1` to the corresponding solver.
 *
 * \param[in] func
 *    The name of the solver for which the size of the workspace should be queried. The argument is not
 *    case sensitive.
 *
 * \param[in] M
 *    The number of rows of the solution, the right hand side, and the order of the matrices from
 *    the left side. If the equation only depends on one dimension argument, this is given here.
 *
 * \param[in] N
 *    The number of columns of the solution, the right hand side, and the order of the matrices from
 *    the right side. If the equation only depends on one dimension argument, this parameter is not used.
 *
 * \return The memory requirements of the routine counted in floating point number of the desired precision.
 * That means the return value either needs to be multiplied by `sizeof(float)` or `sizeof(double)`. A negative
 * return value indicates an error. This error code is explained in the description of the corresponding routine
 * for which the memory query was performed.
 */
ssize_t mepack_memory(const char *func, int M, int N)
{
    Int info = -1;
    ssize_t return_value = -1;

    char *function_name_internal = NULL;

    if (strncasecmp(func, "DLA", 3) == 0 || strncasecmp(func, "SLA", 3) == 0)
    {
        function_name_internal = strdup(func);
    }
    else
    {
        if ( strcasecmp("mepack_single_trlyap_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRLYAP_DAG");
        } else if ( strcasecmp("mepack_double_trlyap_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRLYAP_DAG");
        } else if ( strcasecmp("mepack_single_trlyap_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRLYAP_L3");
        } else if ( strcasecmp("mepack_double_trlyap_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRLYAP_L3");
        } else if ( strcasecmp("mepack_single_trlyap_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRLYAP_L3_2S");
        } else if ( strcasecmp("mepack_double_trlyap_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRLYAP_L3_2S");
        } else if ( strcasecmp("mepack_single_trlyap_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRLYAP_L2");
        } else if ( strcasecmp("mepack_double_trlyap_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRLYAP_L2");
        } else if ( strcasecmp("mepack_single_trlyap_level2_opt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRLYAP_L2_OPT");
        } else if ( strcasecmp("mepack_double_trlyap_level2_opt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRLYAP_L2_OPT");
        } else if ( strcasecmp("mepack_single_trlyap_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRLYAP_RECURSIVE");
        } else if ( strcasecmp("mepack_double_trlyap_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRLYAP_RECURSIVE");
        } else if ( strcasecmp("mepack_single_trstein_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSTEIN_DAG");
        } else if ( strcasecmp("mepack_double_trstein_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSTEIN_DAG");
        } else if ( strcasecmp("mepack_single_trstein_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSTEIN_L3");
        } else if ( strcasecmp("mepack_double_trstein_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSTEIN_L3");
        } else if ( strcasecmp("mepack_single_trstein_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSTEIN_L2");
        } else if ( strcasecmp("mepack_double_trstein_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSTEIN_L2");
        } else if ( strcasecmp("mepack_single_trstein_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSTEIN_L3_2S");
        } else if ( strcasecmp("mepack_double_trstein_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSTEIN_L3_2S");
        } else if ( strcasecmp("mepack_single_trstein_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSTEIN_RECURSIVE");
        } else if ( strcasecmp("mepack_double_trstein_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSTEIN_RECURSIVE");
        } else if ( strcasecmp("mepack_double_tglyap_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGLYAP_DAG");
        } else if ( strcasecmp("mepack_single_tglyap_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGLYAP_DAG");
        } else if ( strcasecmp("mepack_double_tglyap_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGLYAP_L3");
        } else if ( strcasecmp("mepack_single_tglyap_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGLYAP_L3");
        } else if ( strcasecmp("mepack_double_tglyap_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGLYAP_L3_2S");
        } else if ( strcasecmp("mepack_single_tglyap_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGLYAP_L3_2S");
        } else if ( strcasecmp("mepack_double_tglyap_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGLYAP_L2");
        } else if ( strcasecmp("mepack_single_tglyap_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGLYAP_L2");
        } else if ( strcasecmp("mepack_double_tglyap_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGLYAP_L2_UNOPT");
        } else if ( strcasecmp("mepack_single_tglyap_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGLYAP_L2_UNOPT");
        } else if ( strcasecmp("mepack_double_tglyap_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGLYAP_RECURSIVE");
        } else if ( strcasecmp("mepack_single_tglyap_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGLYAP_RECURSIVE");
        } else if ( strcasecmp("mepack_double_tgstein_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSTEIN_DAG");
        } else if ( strcasecmp("mepack_single_tgstein_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSTEIN_DAG");
        } else if ( strcasecmp("mepack_double_tgstein_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSTEIN_L3");
        } else if ( strcasecmp("mepack_single_tgstein_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSTEIN_L3");
        } else if ( strcasecmp("mepack_double_tgstein_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSTEIN_L3_2S");
        } else if ( strcasecmp("mepack_single_tgstein_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSTEIN_L3_2S");
        } else if ( strcasecmp("mepack_double_tgstein_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSTEIN_L2");
        } else if ( strcasecmp("mepack_single_tgstein_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSTEIN_L2");
        } else if ( strcasecmp("mepack_double_tgstein_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSTEIN_RECURSIVE");
        } else if ( strcasecmp("mepack_single_tgstein_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSTEIN_RECURSIVE");
        } else if ( strcasecmp("mepack_double_trsylv_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_DAG");
        } else if ( strcasecmp("mepack_single_trsylv_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_DAG");
        } else if ( strcasecmp("mepack_double_trsylv_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L3");
        } else if ( strcasecmp("mepack_single_trsylv_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L3");
        } else if ( strcasecmp("mepack_double_trsylv_level3_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L3_UNOPT");
        } else if ( strcasecmp("mepack_single_trsylv_level3_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L3_UNOPT");
        } else if ( strcasecmp("mepack_double_trsylv_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L3_2S");
        } else if ( strcasecmp("mepack_single_trsylv_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L3_2S");
        } else if ( strcasecmp("mepack_double_trsylv_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_RECURSIVE");
        } else if ( strcasecmp("mepack_single_trsylv_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_RECURSIVE");
        } else if ( strcasecmp("mepack_double_trsylv_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2");
        } else if ( strcasecmp("mepack_single_trsylv_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2");
        } else if ( strcasecmp("mepack_double_trsylv_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_UNOPT");
        } else if ( strcasecmp("mepack_single_trsylv_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_UNOPT");
        } else if ( strcasecmp("mepack_double_trsylv_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_REORDER");
        } else if ( strcasecmp("mepack_single_trsylv_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_REORDER");
        } else if ( strcasecmp("mepack_double_trsylv_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_single_trsylv_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_double_trsylv_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_single_trsylv_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_double_trsylv_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_single_trsylv_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_double_trsylv_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_single_trsylv_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_double_trsylv_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_single_trsylv_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_double_trsylv2_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_DAG");
        } else if ( strcasecmp("mepack_single_trsylv2_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_DAG");
        } else if ( strcasecmp("mepack_double_trsylv2_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L3");
        } else if ( strcasecmp("mepack_single_trsylv2_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L3");
        } else if ( strcasecmp("mepack_double_trsylv2_level3_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L3_UNOPT");
        } else if ( strcasecmp("mepack_single_trsylv2_level3_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L3_UNOPT");
        } else if ( strcasecmp("mepack_double_trsylv2_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L3_2S");
        } else if ( strcasecmp("mepack_single_trsylv2_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L3_2S");
        } else if ( strcasecmp("mepack_double_trsylv2_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_RECURSIVE");
        } else if ( strcasecmp("mepack_single_trsylv2_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_RECURSIVE");
        } else if ( strcasecmp("mepack_double_trsylv2_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2");
        } else if ( strcasecmp("mepack_single_trsylv2_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_UNOPT");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_UNOPT");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_REORDER");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_REORDER");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_double_trsylv2_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("DLA_TRSYLV2_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_single_trsylv2_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("SLA_TRSYLV2_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_double_tgsylv_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_DAG");
        } else if ( strcasecmp("mepack_single_tgsylv_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_DAG");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_REORDER");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_REORDER");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_rowwise", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_ROWWISE");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_rowwise", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_ROWWISE");
        } else if ( strcasecmp("mepack_double_tgsylv_level2_colwise", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2_COLWISE");
        } else if ( strcasecmp("mepack_single_tgsylv_level2_colwise", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2_COLWISE");
        } else if ( strcasecmp("mepack_double_tgsylv_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L2");
        } else if ( strcasecmp("mepack_single_tgsylv_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L2");
        } else if ( strcasecmp("mepack_double_tgsylv_gardiner_laub", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_GARDINER_LAUB");
        } else if ( strcasecmp("mepack_single_tgsylv_gardiner_laub", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_GARDINER_LAUB");
        } else if ( strcasecmp("mepack_double_tgsylv_level3_rowwise", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L3_ROWWISE");
        } else if ( strcasecmp("mepack_single_tgsylv_level3_rowwise", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L3_ROWWISE");
        } else if ( strcasecmp("mepack_double_tgsylv_level3_colwise", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L3_COLWISE");
        } else if ( strcasecmp("mepack_single_tgsylv_level3_colwise", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L3_COLWISE");
        } else if ( strcasecmp("mepack_double_tgsylv_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L3");
        } else if ( strcasecmp("mepack_single_tgsylv_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L3");
        } else if ( strcasecmp("mepack_double_tgsylv_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_L3_2S");
        } else if ( strcasecmp("mepack_single_tgsylv_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_L3_2S");
        } else if ( strcasecmp("mepack_double_tgsylv_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGSYLV_RECURSIVE");
        } else if ( strcasecmp("mepack_single_tgsylv_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGSYLV_RECURSIVE");
        } else if ( strcasecmp("mepack_double_tgcsylv_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L3");
        } else if ( strcasecmp("mepack_double_tgcsylv_level3_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L3_UNOPT");
        } else if ( strcasecmp("mepack_double_tgcsylv_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L3_2S");
        } else if ( strcasecmp("mepack_double_tgcsylv_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DAG");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_UNOPT");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_REORDER");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_double_tgcsylv_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_double_tgcsylv_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_RECURSIVE");
        } else if ( strcasecmp("mepack_single_tgcsylv_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L3");
        } else if ( strcasecmp("mepack_single_tgcsylv_level3_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L3_UNOPT");
        } else if ( strcasecmp("mepack_single_tgcsylv_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L3_2S");
        } else if ( strcasecmp("mepack_single_tgcsylv_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DAG");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_unopt", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_UNOPT");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_reorder", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_REORDER");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_single_tgcsylv_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_single_tgcsylv_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_RECURSIVE");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level3", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L3");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L3_2S");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_dag", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_DAG");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level2", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L2");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_double_tgcsylv_dual_recursive", func ) == 0 ) {
            function_name_internal = strdup("DLA_TGCSYLV_DUAL_RECURSIVE");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level3", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L3");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level3_2stage", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L3_2S");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_dag", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_DAG");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level2", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L2");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level2_local_copy", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level2_local_copy_32", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level2_local_copy_64", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level2_local_copy_96", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_level2_local_copy_128", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128");
        } else if ( strcasecmp("mepack_single_tgcsylv_dual_recursive", func ) == 0 ) {
            function_name_internal = strdup("SLA_TGCSYLV_DUAL_RECURSIVE");
        }

    }
    /* Check for Local Copies */
    if (strstr(function_name_internal, "LOCAL_COPY_32") != NULL
            && MAX(M,N) > 32) {
        return_value = -1;
        goto end;
    }

    if (strstr(function_name_internal, "LOCAL_COPY_64") != NULL
            && MAX(M,N) > 64) {
        return_value = -1;
        goto end;
    }

    if (strstr(function_name_internal, "LOCAL_COPY_96") != NULL
            && MAX(M,N) > 96) {
        return_value = -1;
        goto end;
    }
    if (strstr(function_name_internal, "LOCAL_COPY_128") != NULL
            && MAX(M,N) > 128) {
        return_value = -1;
        goto end;
    }
    if (strstr(function_name_internal, "LOCAL_COPY") != NULL
            && MAX(M,N) > 128) {
        return_value = -1;
        goto end;
    }


    /*
     * TRLYAP (DOUBLE)
     */
    if ((strcasecmp("DLA_TRLYAP_L3", function_name_internal) == 0)) {
        FC_GLOBAL_(dla_trlyap_l3,DLA_TRLYAP_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TRLYAP_L3_2S", function_name_internal) == 0)) {
        FC_GLOBAL_(dla_trlyap_l3_2s,DLA_TRLYAP_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TRLYAP_DAG", function_name_internal) == 0)){
        FC_GLOBAL_(dla_trlyap_dag,DLA_TRLYAP_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TRLYAP_L2_OPT", function_name_internal) == 0)) {
        FC_GLOBAL_(dla_trlyap_l2_opt,DLA_TRLYAP_L2_OPT)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TRLYAP_L2", function_name_internal) == 0) ) {
        FC_GLOBAL_(dla_trlyap_l2,DLA_TRLYAP_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TRLYAP_RECURSIVE", function_name_internal) == 0))  {
        FC_GLOBAL_(dla_trlyap_recursive,DLA_TRLYAP_RECURSIVE)("N", &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TRLYAP (SINGLE)
     */
    if ((strcasecmp("SLA_TRLYAP_L3", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_trlyap_l3,SLA_TRLYAP_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TRLYAP_L3_2S", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_trlyap_l3_2s,SLA_TRLYAP_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TRLYAP_DAG", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_trlyap_dag,SLA_TRLYAP_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TRLYAP_L2_OPT", function_name_internal) == 0))  {
        FC_GLOBAL_(sla_trlyap_l2_opt,SLA_TRLYAP_L2_OPT)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TRLYAP_L2", function_name_internal) == 0) ) {
        FC_GLOBAL_(sla_trlyap_l2,SLA_TRLYAP_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TRLYAP_RECURSIVE", function_name_internal) == 0)){
        FC_GLOBAL_(sla_trlyap_recursive,SLA_TRLYAP_RECURSIVE)("N", &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }



    /*
     * TRSTEIN (DOUBLE)
     */
    if ((strcasecmp("DLA_TRSTEIN_L3", function_name_internal) == 0)  ) {
        FC_GLOBAL_(dla_trstein_l3,DLA_TRSTEIN_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TRSTEIN_L3_2S", function_name_internal) == 0) ) {
        FC_GLOBAL_(dla_trstein_l3_2s,DLA_TRSTEIN_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TRSTEIN_DAG", function_name_internal) == 0)  ) {
        FC_GLOBAL_(dla_trstein_dag,DLA_TRSTEIN_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TRSTEIN_L2", function_name_internal) == 0)  ) {
        FC_GLOBAL_(dla_trstein_l2,DLA_TRSTEIN_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TRSTEIN_RECURSIVE", function_name_internal) == 0)  ) {
        FC_GLOBAL_(dla_trstein_recursive,DLA_TRSTEIN_RECURSIVE)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TRSTEIN (SINGLE)
     */
    if ((strcasecmp("SLA_TRSTEIN_L3", function_name_internal) == 0)  ) {
        FC_GLOBAL_(sla_trstein_l3,SLA_TRSTEIN_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TRSTEIN_L3_2S", function_name_internal) == 0) ) {
        FC_GLOBAL_(sla_trstein_l3_2s,SLA_TRSTEIN_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TRSTEIN_DAG", function_name_internal) == 0)  ) {
        FC_GLOBAL_(sla_trstein_dag,SLA_TRSTEIN_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TRSTEIN_L2", function_name_internal) == 0) ) {
        FC_GLOBAL_(sla_trstein_l2,SLA_TRSTEIN_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TRSTEIN_RECURSIVE", function_name_internal) == 0)  ) {
        FC_GLOBAL_(sla_trstein_recursive,SLA_TRSTEIN_RECURSIVE)("N",  &(Int){M}, NULL /*A*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TRSYLV (DOUBLE)
     */
    if ( strcasecmp("DLA_TRSYLV_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l3,DLA_TRSYLV_L3)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L3_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l3_unopt,DLA_TRSYLV_L3_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_dag,DLA_TRSYLV_DAG)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l3_2s,DLA_TRSYLV_L3_2S)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2,DLA_TRSYLV_L2)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_unopt,DLA_TRSYLV_L2_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_REORDER", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_reorder,DLA_TRSYLV_L2_REORDER)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_local_copy,DLA_TRSYLV_L2_LOCAL_COPY)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_32", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_local_copy_32,DLA_TRSYLV_L2_LOCAL_COPY_32)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_64", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_local_copy_64,DLA_TRSYLV_L2_LOCAL_COPY_64)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_96", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_local_copy_96,DLA_TRSYLV_L2_LOCAL_COPY_96)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_128", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_l2_local_copy_128,DLA_TRSYLV_L2_LOCAL_COPY_128)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv_recursive,DLA_TRSYLV_RECURSIVE)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TRSYLV (Single)
     */
    if ( strcasecmp("SLA_TRSYLV_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l3,SLA_TRSYLV_L3)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L3_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l3_unopt,SLA_TRSYLV_L3_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_dag,SLA_TRSYLV_DAG)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l3_2s,SLA_TRSYLV_L3_2S)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2,SLA_TRSYLV_L2)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_unopt,SLA_TRSYLV_L2_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_REORDER", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_reorder,SLA_TRSYLV_L2_REORDER)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_local_copy,SLA_TRSYLV_L2_LOCAL_COPY)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_32", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_local_copy_32,SLA_TRSYLV_L2_LOCAL_COPY_32)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_64", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_local_copy_64,SLA_TRSYLV_L2_LOCAL_COPY_64)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_96", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_local_copy_96,SLA_TRSYLV_L2_LOCAL_COPY_96)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_128", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_l2_local_copy_128,SLA_TRSYLV_L2_LOCAL_COPY_128)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv_recursive,SLA_TRSYLV_RECURSIVE)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TRSYLV2 (DOUBLE)
     */
    if ( strcasecmp("DLA_TRSYLV2_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l3,DLA_TRSYLV2_L3)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L3_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l3_unopt,DLA_TRSYLV2_L3_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_dag,DLA_TRSYLV2_DAG)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l3_2s,DLA_TRSYLV2_L3_2S)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2,DLA_TRSYLV2_L2)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_unopt,DLA_TRSYLV2_L2_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_REORDER", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_reorder,DLA_TRSYLV2_L2_REORDER)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_local_copy,DLA_TRSYLV2_L2_LOCAL_COPY)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_32", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_local_copy_32,DLA_TRSYLV2_L2_LOCAL_COPY_32)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_64", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_local_copy_64,DLA_TRSYLV2_L2_LOCAL_COPY_64)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_96", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_local_copy_96,DLA_TRSYLV2_L2_LOCAL_COPY_96)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_128", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_l2_local_copy_128,DLA_TRSYLV2_L2_LOCAL_COPY_128)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TRSYLV2_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_trsylv2_recursive,DLA_TRSYLV2_RECURSIVE)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TRSYLV2 (Single)
     */
    if ( strcasecmp("SLA_TRSYLV2_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l3,SLA_TRSYLV2_L3)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L3_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l3_unopt,SLA_TRSYLV2_L3_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_dag,SLA_TRSYLV2_DAG)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l3_2s,SLA_TRSYLV2_L3_2S)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2,SLA_TRSYLV2_L2)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_unopt,SLA_TRSYLV2_L2_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_REORDER", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_reorder,SLA_TRSYLV2_L2_REORDER)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_local_copy,SLA_TRSYLV2_L2_LOCAL_COPY)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_32", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_local_copy_32,SLA_TRSYLV2_L2_LOCAL_COPY_32)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_64", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_local_copy_64,SLA_TRSYLV2_L2_LOCAL_COPY_64)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_96", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_local_copy_96,SLA_TRSYLV2_L2_LOCAL_COPY_96)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_128", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_l2_local_copy_128,SLA_TRSYLV2_L2_LOCAL_COPY_128)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TRSYLV2_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_trsylv2_recursive,SLA_TRSYLV2_RECURSIVE)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }


    /*
     * TGLYAP (DOUBLE)
     */
    if ( strcasecmp("DLA_TGLYAP_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l3,DLA_TGLYAP_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGLYAP_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l3_2s,DLA_TGLYAP_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGLYAP_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tglyap_dag,DLA_TGLYAP_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGLYAP_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l2,DLA_TGLYAP_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGLYAP_L2_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l2_unopt,DLA_TGLYAP_L2_UNOPT)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGLYAP_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tglyap_recursive,DLA_TGLYAP_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }




    /*
     * TGLYAP (SINGLE)
     */
    if ( strcasecmp("SLA_TGLYAP_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l3,SLA_TGLYAP_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGLYAP_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l3_2s,SLA_TGLYAP_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGLYAP_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tglyap_dag,SLA_TGLYAP_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGLYAP_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l2,SLA_TGLYAP_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGLYAP_L2_UNOPT", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l2_unopt,SLA_TGLYAP_L2_UNOPT)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGLYAP_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tglyap_recursive,SLA_TGLYAP_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TGSTEIN (DOUBLE)
     */
    if ( strcasecmp("DLA_TGSTEIN_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tgstein_l3,DLA_TGSTEIN_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGSTEIN_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tgstein_l3_2s,DLA_TGSTEIN_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGSTEIN_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tgstein_dag,DLA_TGSTEIN_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGSTEIN_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tgstein_l2,DLA_TGSTEIN_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGSTEIN_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_tgstein_recursive,DLA_TGSTEIN_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }



    /*
     * TGSTEIN(SINGLE)
     */
    if ( strcasecmp("SLA_TGSTEIN_L3", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tgstein_l3,SLA_TGSTEIN_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGSTEIN_L3_2S", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tgstein_l3_2s,SLA_TGSTEIN_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGSTEIN_DAG", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tgstein_dag,SLA_TGSTEIN_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGSTEIN_L2", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tgstein_l2,SLA_TGSTEIN_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGSTEIN_RECURSIVE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_tgstein_recursive,SLA_TGSTEIN_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/ , &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return_value = (ssize_t) info;
        goto end;
    }





    /*
     * TGSYLV (DOUBLE)
     */
    if ((strcasecmp("DLA_TGSYLV_L3", function_name_internal ) == 0 ) ) {
        FC_GLOBAL_(dla_tgsylv_l3,DLA_TGSYLV_L3)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L3_ROWWISE", function_name_internal ) == 0 )) {
        FC_GLOBAL_(dla_tgsylv_l3_rowwise,DLA_TGSYLV_L3_ROWWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L3_COLWISE", function_name_internal ) == 0 )) {
        FC_GLOBAL_(dla_tgsylv_l3_colwise,DLA_TGSYLV_L3_COLWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L3_2S", function_name_internal ) == 0 )) {
        FC_GLOBAL_(dla_tgsylv_l3_2s,DLA_TGSYLV_L3_2S)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_DAG", function_name_internal ) == 0 ) ) {
        FC_GLOBAL_(dla_tgsylv_dag,DLA_TGSYLV_DAG)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_32", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_32,DLA_TGSYLV_L2_LOCAL_COPY_32)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_64", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_64,DLA_TGSYLV_L2_LOCAL_COPY_64)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_96", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_96,DLA_TGSYLV_L2_LOCAL_COPY_96)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_128", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_128,DLA_TGSYLV_L2_LOCAL_COPY_128)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy,DLA_TGSYLV_L2_LOCAL_COPY)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_REORDER", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_reorder,DLA_TGSYLV_L2_REORDER)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_L2", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2,DLA_TGSYLV_L2)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_ROWWISE", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_rowwise,DLA_TGSYLV_L2_ROWWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_COLWISE", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_colwise,DLA_TGSYLV_L2_COLWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_UNOPT", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_unopt,DLA_TGSYLV_L2_UNOPT)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_RECURSIVE", function_name_internal ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_recursive,DLA_TGSYLV_RECURSIVE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("DLA_TGSYLV_GARDINER_LAUB", function_name_internal ) == 0 ) || (strcasecmp("SLA_TGSYLV_GARDINER_LAUB", function_name_internal) == 0)) {
        FC_GLOBAL_(dla_tgsylv_gardiner_laub,DLA_TGSYLV_GARDINER_LAUB)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TGSYLV (SINGLE)
     */
    if ((strcasecmp("SLA_TGSYLV_L3", function_name_internal ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3,SLA_TGSYLV_L3)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L3_ROWWISE", function_name_internal ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3_ROWISE", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3_rowwise,SLA_TGSYLV_L3_ROWWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L3_COLWISE", function_name_internal ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3_COLWISE", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3_colwise,SLA_TGSYLV_L3_COLWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L3_2S", function_name_internal ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3_2S", function_name_internal) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3_2s,SLA_TGSYLV_L3_2S)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_DAG", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_dag,SLA_TGSYLV_DAG)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_32", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_32,SLA_TGSYLV_L2_LOCAL_COPY_32)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_64", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_64,SLA_TGSYLV_L2_LOCAL_COPY_64)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_96", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_96,SLA_TGSYLV_L2_LOCAL_COPY_96)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_128", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_128,SLA_TGSYLV_L2_LOCAL_COPY_128)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy,SLA_TGSYLV_L2_LOCAL_COPY)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_REORDER", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_reorder,SLA_TGSYLV_L2_REORDER)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_L2", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2,SLA_TGSYLV_L2)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_ROWWISE", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_rowwise,SLA_TGSYLV_L2_ROWWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_COLWISE", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_colwise,SLA_TGSYLV_L2_COLWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_UNOPT", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_unopt,SLA_TGSYLV_L2_UNOPT)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_RECURSIVE", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_recursive,SLA_TGSYLV_RECURSIVE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    if ((strcasecmp("SLA_TGSYLV_GARDINER_LAUB", function_name_internal ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_gardiner_laub,SLA_TGSYLV_GARDINER_LAUB)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TGCSYLV(DOUBLE)
     */
    if ( strcasecmp("DLA_TGCSYLV_L2",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2,DLA_TGCSYLV_L2)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_UNOPT",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_unopt,DLA_TGCSYLV_L2_UNOPT)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_REORDER",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_reorder,DLA_TGCSYLV_L2_REORDER)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy,DLA_TGCSYLV_L2_LOCAL_COPY)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_32",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_32,DLA_TGCSYLV_L2_LOCAL_COPY_32)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_64",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_64,DLA_TGCSYLV_L2_LOCAL_COPY_64)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_96",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_96,DLA_TGCSYLV_L2_LOCAL_COPY_96)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_128",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_128,DLA_TGCSYLV_L2_LOCAL_COPY_128)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L3",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l3,DLA_TGCSYLV_L3)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L3_UNOPT",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l3_unopt,DLA_TGCSYLV_L3_UNOPT)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_L3_2S",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l3_2s,DLA_TGCSYLV_L3_2S)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DAG",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dag,DLA_TGCSYLV_DAG)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_RECURSIVE",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_recursive,DLA_TGCSYLV_RECURSIVE)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TGCSYLV_DUAL_DUAL(DOUBLE)
     */
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2,DLA_TGCSYLV_DUAL_L2)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_32,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_64,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_96,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_128,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L3",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l3,DLA_TGCSYLV_DUAL_L3)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L3_2S",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l3_2s,DLA_TGCSYLV_DUAL_L3_2S)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_DAG",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_dag,DLA_TGCSYLV_DUAL_DAG)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_RECURSIVE",function_name_internal) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_recursive,DLA_TGCSYLV_DUAL_RECURSIVE)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TGCSYLV(SINGLE)
     */
    if ( strcasecmp("SLA_TGCSYLV_L2",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2,SLA_TGCSYLV_L2)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_UNOPT",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_unopt,SLA_TGCSYLV_L2_UNOPT)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_REORDER",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_reorder,SLA_TGCSYLV_L2_REORDER)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy,SLA_TGCSYLV_L2_LOCAL_COPY)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_32",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_32,SLA_TGCSYLV_L2_LOCAL_COPY_32)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_64",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_64,SLA_TGCSYLV_L2_LOCAL_COPY_64)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_96",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_96,SLA_TGCSYLV_L2_LOCAL_COPY_96)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_128",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_128,SLA_TGCSYLV_L2_LOCAL_COPY_128)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L3",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l3,SLA_TGCSYLV_L3)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L3_UNOPT",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l3_unopt,SLA_TGCSYLV_L3_UNOPT)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_L3_2S",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l3_2s,SLA_TGCSYLV_L3_2S)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DAG",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dag,SLA_TGCSYLV_DAG)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_RECURSIVE",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_recursive,SLA_TGCSYLV_RECURSIVE)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }

    /*
     * TGCSYLV_DUAL(SINGLE)
     */
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2,SLA_TGCSYLV_DUAL_L2)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_32,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_64,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_96,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_128,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L3",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l3,SLA_TGCSYLV_DUAL_L3)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L3_2S",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l3_2s,SLA_TGCSYLV_DUAL_L3_2S)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_DAG",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_dag,SLA_TGCSYLV_DUAL_DAG)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_RECURSIVE",function_name_internal) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_recursive,SLA_TGCSYLV_DUAL_RECURSIVE)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int) {N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int) {N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return_value = (ssize_t) info;
        goto end;
    }



end:
    free(function_name_internal);
    return return_value;
}


/**
 * @}
 */
